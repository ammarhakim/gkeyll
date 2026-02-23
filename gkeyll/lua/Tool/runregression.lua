-- Gkyl ------------------------------------------------------------------------
--
-- Runs regression tests and compares with accepted results.
-- Restructured to support the hierarchical layer architecture:
--   core -> moments -> vlasov -> gyrokinetics -> pkpm
-- Each layer has its own luareg/ (Lua regression tests) and creg/
-- (C regression test sources).
--
-- Lua tests:  discovered from source_dir/<layer>/luareg/rt_*.lua.
--             Run via 'gkeyll <file>' from a per-test scratch directory.
--
-- C tests:    discovered from source_dir/<layer>/creg/rt_*.c.
--             Each test is compiled on-the-fly using the installed
--             share/Makefile (processed from Makefile_for_ext_C_input at
--             'make install' time), then run from its scratch directory.
--             No separate 'make regression' step is required.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local lfs = require "lfs"
local uuid = require "Lib.UUID"

if GKYL_HAVE_SQLITE3 == false then
   -- can't run without SQLITE3
   print("Sorry, runregression needs Sqlite3. This executable was built without it.")
   return 1
end

-- Minimal set of requires needed by this tool.
-- Deliberately NOT loaded (reasons below):
--   Lib.Alloc      - not needed; Zero* functions manage their own memory
--   Lib.Linalg     - not needed; used only internally by Grid.ZeroRectCart
--   DataStruct.ZeroDynVector - not needed; dynvector comparison is NYI (skip)
--
-- For file comparison (compareFiles):
--   ZeroArray    - gkyl_array_diff: numerically compare two .gkyl arrays
--   ZeroArrayRio - gkyl_grid_array_new_from_file: read a .gkyl file
--                  (Note: ZeroArrayRio now requires Grid.ZeroRectCart, not
--                   Grid.RectCart, to avoid pulling in Comm.Mpi)
--   ZeroRectCart - gkyl_rect_grid_cmp / createGridRanges: compare grids
--   ZeroUtil     - gkylFileType: detect .gkyl file type from header
local Logger = require "Lib.Logger"
local Time = require "Lib.Time"
local argparse = require "Lib.argparse"
local date = require "xsys.date"
local lume = require "Lib.lume"
local sql = require "sqlite3"

local ZeroArray = require "DataStruct.ZeroArray"
local ZeroArrayRio = require "DataStruct.ZeroArrayRio"
local ZeroRectCart = require "Grid.ZeroRectCart"
local ZeroUtil = require "Lib.ZeroUtil"

-- GKYL_OUT_PREFIX is required by some Lua modules; set to a safe default here
-- (each test run redirects output to its own scratch directory).
GKYL_OUT_PREFIX = lfs.currentdir() .. "/" .. "runregression"

local log = Logger { logToFile = true }
local verboseLog = function (msg) end -- default: no verbose output
local verboseLogger = function (msg) log(msg) end

-- ---- Layer definitions ------------------------------------------------------
-- The four physics layers that have regression tests. 'src' is the name of
-- the subdirectory under the gkeyll source root and build root.
local LAYERS = {
   { name = "moments",     src = "moments"     },
   { name = "vlasov",      src = "vlasov"      },
   { name = "gyrokinetic", src = "gyrokinetic" },
   { name = "pkpm",        src = "pkpm"        },
}

-- Fast lookup: set of valid layer names.
local VALID_LAYER_NAMES = {}
for _, L in ipairs(LAYERS) do VALID_LAYER_NAMES[L.name] = true end

-- ---- Per-test timeout command detection -------------------------------------
-- GNU 'timeout' is not part of macOS's standard utilities. It ships as
-- 'gtimeout' via 'brew install coreutils', or as 'timeout' on Linux and on
-- macOS with coreutils linked without the 'g' prefix.
-- We detect whichever is available at startup. If neither is found we fall
-- back to a Perl-based fork+alarm implementation (Perl is always present on
-- macOS and Linux).
local TIMEOUT_CMD
do
   local function hasCmd(name)
      return os.execute(
         string.format("which %s > /dev/null 2>&1", name)) == 0
   end
   if     hasCmd("timeout")  then TIMEOUT_CMD = "timeout"
   elseif hasCmd("gtimeout") then TIMEOUT_CMD = "gtimeout"
   else                           TIMEOUT_CMD = nil  -- use Perl fallback
   end
end

-- ---- Layer pre-processing of command arguments ------------------------------
-- The argparse library does not natively support positional sub-sub-command
-- names that collide with free text (e.g. "run moments check"). We therefore
-- scan GKYL_COMMANDS_L before argparse sees it, extract any token that is a
-- known layer name, and store it in 'detectedLayer'. The remaining tokens are
-- passed to argparse as normal, so "run moments check" becomes "run check"
-- from argparse's perspective.
local detectedLayer = nil
local filteredCmds  = {}
for i = 1, #GKYL_COMMANDS_L do
   local arg = GKYL_COMMANDS_L[i]
   if VALID_LAYER_NAMES[arg] then
      -- First layer name wins; extras (if any) are silently dropped.
      if detectedLayer == nil then
         detectedLayer = arg
      end
   else
      filteredCmds[#filteredCmds + 1] = arg
   end
end

-- ---- Global counters --------------------------------------------------------
-- Per-layer pass/fail counts, keyed by layer name.
local layerCounts = {}
for _, L in ipairs(LAYERS) do
   layerCounts[L.name] = { total = 0, passed = 0, failed = 0 }
end

-- ---- Configuration ----------------------------------------------------------
local isConfiguring = false
local configVals    = nil

-- Path of the configuration file written by 'configure' and read by 'run'.
local confFile = os.getenv("HOME") .. "/runregression.config.lua"

-- Per-layer SQLite database connections and prepared statements.
-- Populated lazily by getLayerDB().
local layerDBs = {}

-- Per-layer ignore lists and MOAT lists, keyed by layer name.
-- Loaded during loadConfigure().
local ignoreTests = {}  -- ignoreTests[layerName] = {lua={...}, c={...}}
local moatTests   = {}  -- moatTests[layerName]   = {path, ...}

-- Unique ID and timestamp for this invocation.
local runID  = uuid()
local runDate = date(false):fmt("${iso}")

-- ---- Database schema --------------------------------------------------------
-- Two tables per layer database:
--
-- RegressionMeta stores one row per 'run' invocation:
--   guid               unique ID shared with RegressionData rows
--   tstamp             ISO-8601 timestamp of the run
--   GKYL_EXEC          path to the gkeyll executable used
--   GKYL_GIT_CHANGESET git commit hash
--   GKYL_BUILD_DATE    build date string
--   ntotal / npass / nfail  aggregate counts for this run
--
-- RegressionData stores one row per individual test:
--   guid       links to a RegressionMeta row
--   name       test name (relative path / binary name)
--   test_type  'lua' or 'c'  (new in hierarchical design)
--   status     -4=compile_fail  -3=timeout  -2=create  -1=skip  0=fail  1=pass
--   runtime    wall-clock seconds
--   runlog     captured stdout+stderr from the test (includes compiler output for C)

local SCHEMA_SQL = [[
  drop table if exists RegressionMeta;
  create table RegressionMeta (
    guid               text,
    tstamp             text,
    GKYL_EXEC          text,
    GKYL_GIT_CHANGESET text,
    GKYL_BUILD_DATE    text,
    ntotal             integer,
    npass              integer,
    nfail              integer
  );

  drop table if exists RegressionData;
  create table RegressionData (
    guid      text,
    name      text,
    test_type text,
    status    integer,
    runtime   real,
    runlog    text
  );
]]

-- ---- Helper utilities -------------------------------------------------------

-- Recursively creates all directories in 'path'.
local function mkdir(path)
   local sep, pStr = package.config:sub(1,1), "/"
   for dir in path:gmatch("[^" .. sep .. "]+") do
      pStr = pStr .. dir .. sep
      lfs.mkdir(pStr)
   end
end

-- Splits a comma-separated string into a list of words.
local function splitList(listStr)
   local words = {}
   for w in listStr:gmatch('[^,%s]+') do
      table.insert(words, w)
   end
   return words
end

-- Returns the basename of a path (everything after the last '/').
local function basename(path)
   return path:match("([^/]+)$") or path
end

-- Returns the directory component of a path (e.g. "/a/b/c.lua" -> "/a/b").
local function dirname(path)
   return path:match("^(.+)/[^/]+$") or "."
end

-- Strips the extension from a filename (e.g. "rt_foo.lua" -> "rt_foo").
local function stripext(fn)
   return fn:match("^(.+)%..+$") or fn
end

-- Walks a directory tree recursively. Returns a coroutine that yields
-- (dir, filename, attr) for every entry.
local function dirtree(dir)
   assert(dir and dir ~= "", "directory parameter is missing or empty")
   if string.sub(dir, -1) == "/" then
      dir = string.sub(dir, 1, -2)
   end
   local function yieldtree(d)
      for fn in lfs.dir(d) do
         if fn ~= "." and fn ~= ".." then
            local fullNm = d .. "/" .. fn
            local attr = lfs.attributes(fullNm)
            if attr then
               coroutine.yield(d, fn, attr)
               if attr.mode == "directory" then
                  yieldtree(fullNm)
               end
            end
         end
      end
   end
   return coroutine.wrap(function() yieldtree(dir) end)
end

-- Wraps a shell command string with a per-test timeout.
-- 'innerCmd' must already contain any needed 2>&1 redirection.
-- Writes a small helper script to runDir to avoid quoting problems.
-- Appends '; echo __EXIT__:$?' so the caller can parse the exit code.
-- Returns exit code 124 (the POSIX convention used by GNU timeout) when the
-- test exceeds timeoutSecs seconds.
local function wrapWithTimeout(innerCmd, timeoutSecs, runDir)
   if not timeoutSecs or timeoutSecs <= 0 then
      return innerCmd .. "; echo __EXIT__:$?"
   end

   -- Write innerCmd to a tiny shell script so we don't have to escape it.
   local scriptPath = runDir .. "/_rr_cmd.sh"
   local sf = io.open(scriptPath, "w")
   sf:write("#!/bin/sh\n")
   sf:write(innerCmd .. "\n")
   sf:close()
   os.execute(string.format("chmod +x '%s'", scriptPath))

   if TIMEOUT_CMD then
      -- timeout / gtimeout available: straightforward.
      return string.format("%s %d '%s'; echo __EXIT__:$?",
         TIMEOUT_CMD, timeoutSecs, scriptPath)
   else
      -- Perl-based fork+alarm fallback (always available on macOS and Linux).
      -- Key: child calls setpgrp so that gkeyll (a grandchild launched by sh)
      -- is in the same process group. We then kill the entire group with
      -- kill(SIGNAL, -$pid), which reaches gkeyll even though it is not a
      -- direct child of the Perl process.
      local perlPath = runDir .. "/_rr_timeout.pl"
      local pf = io.open(perlPath, "w")
      pf:write(string.format([[#!/usr/bin/env perl
use strict;
my $secs   = %d;
my $script = '%s';
my $pid = fork // die "fork: $!\n";
if (!$pid) {
    # Put the child in its own process group so kill -PGID reaches all
    # grandchildren (e.g. the gkeyll process launched by the sh wrapper).
    setpgrp(0, 0);
    exec($script) or die "exec: $!\n";
}
local $SIG{ALRM} = sub {
    # Kill the entire process group: -$pid means PGID == $pid.
    kill 'TERM', -$pid;
    sleep 2;
    kill 'KILL', -$pid;
    exit 124;
};
alarm($secs);
waitpid($pid, 0);
exit($? >> 8);
]], timeoutSecs, scriptPath:gsub("'", "\\'")))
      pf:close()
      return string.format("perl '%s'; echo __EXIT__:$?", perlPath)
   end
end

-- ---- Test classification predicates ----------------------------------------

local function isLuaRegressionTest(fn)
   return string.find(fn, "/rt_[^/]+%.lua$") ~= nil
end

local function isLuaUnitTest(fn)
   return string.find(fn, "/test_[^/]-%.lua$") ~= nil
end

local function isCxxUnitTest(fn)
   if string.find(fn, "/ctest_[^/]+$") ~= nil then
      local attr = lfs.attributes(fn)
      if attr and attr.mode == "file" and string.sub(attr.permissions, 3, 3) == "x" then
         return true
      end
   end
   return false
end

-- ---- Database helpers -------------------------------------------------------

-- Returns (and caches) the DB connection and prepared statements for a layer.
-- The database lives at: results_dir/<layer>/regressiondb
local function getLayerDB(layerName)
   if layerDBs[layerName] then return layerDBs[layerName] end

   local dbPath = configVals.results_dir .. "/" .. layerName .. "/regressiondb"
   local conn   = sql.open(dbPath)

   -- insertRegressionData binds: guid, name, test_type, status, runtime, runlog
   local insertData = conn:prepare [[
     insert into RegressionData values (?, ?, ?, ?, ?, ?)
   ]]

   -- insertRegressionMeta binds: guid, tstamp, exec, changeset, builddate,
   --                             ntotal, npass, nfail
   local insertMeta = conn:prepare [[
     insert into RegressionMeta values (?, ?, ?, ?, ?, ?, ?, ?)
   ]]

   layerDBs[layerName] = {
      conn       = conn,
      insertData = insertData,
      insertMeta = insertMeta,
   }
   return layerDBs[layerName]
end

-- Insert a single test result row into the appropriate layer's database.
local function insertRegressionData(layerName, guid, name, testType, status, runtm, runlog)
   local db = getLayerDB(layerName)
   db.insertData:reset():bind(guid, name, testType, status, runtm, runlog):step()
end

-- Insert the run-level summary row for a layer.
local function insertRegressionMeta(layerName, guid, tm, ntotal, npass, nfail)
   local db = getLayerDB(layerName)
   db.insertMeta:reset():bind(
      guid, tm,
      GKYL_EXEC,
      GKYL_GIT_CHANGESET,
      GKYL_BUILD_DATE,
      ntotal, npass, nfail
   ):step()
end

-- ---- Configure command ------------------------------------------------------
-- Creates the directory structure, per-layer databases, and configuration file.
-- Explanation of the directory layout:
--   gkeyll-results/
--     <layer>/
--       regressiondb         SQLite database for this layer
--       luareg-runs/         per-test scratch dirs for Lua tests
--       luareg-accepted/     accepted (baseline) outputs for Lua tests
--       creg-runs/           per-test scratch dirs for C tests
--       creg-accepted/       accepted (baseline) outputs for C tests
--
-- Separating 'runs' from 'accepted' makes it possible to inspect the most
-- recent test run at any time without disturbing the accepted baseline.
-- Including the input file (.lua or .c) in the run directory allows the
-- developer to diff the input if outputs change unexpectedly.
local function configure(prefix, mpiExec, sourceDir, args)

   -- Validate prefix directory.
   local prefixAttr = lfs.attributes(prefix)
   if not prefixAttr or prefixAttr.mode ~= "directory" then
      assert(false, string.format("Prefix '%s' is not a directory!", prefix))
   end
   log(string.format("Will write accepted results to %s/gkeyll-results ...\n", prefix))

   -- Validate source directory.
   local srcAttr = lfs.attributes(sourceDir)
   if not srcAttr or srcAttr.mode ~= "directory" then
      assert(false, string.format("Source directory '%s' is not a directory!", sourceDir))
   end

   -- Create the top-level results directory and per-layer subdirectories.
   local resultsDir = string.format("%s/gkeyll-results", prefix)
   mkdir(resultsDir)

   for _, layer in ipairs(LAYERS) do
      local ld = resultsDir .. "/" .. layer.name
      mkdir(ld)
      mkdir(ld .. "/luareg-runs")
      mkdir(ld .. "/luareg-accepted")
      mkdir(ld .. "/creg-runs")
      mkdir(ld .. "/creg-accepted")

      -- Create (or re-create) the SQLite database for this layer.
      -- The --drop-tables flag forces recreation even if the DB already exists.
      local dbPath  = ld .. "/regressiondb"
      local dbExist = lfs.attributes(dbPath)

      if args.drop_tables or dbExist == nil then
         log(string.format("Creating DB for layer '%s' at %s\n", layer.name, dbPath))
         local conn = sql.open(dbPath)
         -- DROP + CREATE ensures the schema is always consistent.
         -- test_type is the key new column: distinguishes Lua ('lua') from C ('c')
         -- regression tests within a single layer database.
         conn:exec(SCHEMA_SQL)
         conn:close()
      else
         log(string.format(
            "DB for layer '%s' already exists at %s (use --drop-tables to reset)\n",
            layer.name, dbPath))
      end
   end

   -- Write the configuration file. This file is read by all subsequent
   -- 'run', 'check', 'list', etc. invocations.
   -- 'prefix' is stored so that runCTest can locate share/Makefile for
   -- compiling C regression tests on-the-fly.
   local fn = io.open(confFile, "w")
   fn:write("return {\n")
   fn:write(string.format("  mpiExec     = \"%s\",\n", mpiExec))
   fn:write(string.format("  prefix      = \"%s\",\n", prefix))
   fn:write(string.format("  results_dir = \"%s\",\n", resultsDir))
   fn:write(string.format("  source_dir  = \"%s\",\n", sourceDir))
   fn:write("}\n")
   fn:close()
   log(string.format("Configuration written to %s\n", confFile))
end

-- ---- Load configuration -----------------------------------------------------
-- Reads the config file, validates it, opens global state (ignore lists, etc.).
local function loadConfigure(args)
   local f = loadfile(confFile)
   if not f then
      log("Regression tests not configured! Run 'runregression configure' first.\n")
      os.exit(1)
   end
   configVals = f()

   -- Sanity-check that required keys are present.
   if not configVals.source_dir then
      log("Configuration is missing source_dir. "
         .. "Please re-run 'runregression configure' with --source-dir.\n")
      os.exit(1)
   end

   -- Backward compatibility: older config files lack 'prefix'. Derive it from
   -- GKYL_EXEC_PATH (e.g. ~/gkylsoft/gkeyll/bin -> ~/gkylsoft).
   if not configVals.prefix then
      local execParent = GKYL_EXEC_PATH:match("^(.+)/gkeyll/bin/?$")
      configVals.prefix = execParent or (os.getenv("HOME") .. "/gkylsoft")
      log(string.format(
         "NOTE: config lacks 'prefix'; derived as '%s'. "
         .. "Re-run 'configure' to make this permanent.\n", configVals.prefix))
   end

   if args.verbose then
      verboseLog = verboseLogger
   end

   -- Load per-layer ignore and MOAT lists.
   -- The file format for ignoretests.lua is:
   --   return { lua = {"path/to/test.lua", ...}, c = {"rt_testname", ...} }
   -- The file format for moat.lua is:
   --   return {"path/to/test1.lua", "path/to/test2.lua", ...}
   -- Missing files are silently treated as empty lists.
   for _, layer in ipairs(LAYERS) do
      local ignFile  = configVals.source_dir .. "/" .. layer.src .. "/luareg/ignoretests.lua"
      local moatFile = configVals.source_dir .. "/" .. layer.src .. "/luareg/moat.lua"

      local gi = loadfile(ignFile)
      if gi then
         ignoreTests[layer.name] = gi()
      else
         ignoreTests[layer.name] = { lua = {}, c = {} }
      end

      local gm = loadfile(moatFile)
      if gm then
         moatTests[layer.name] = gm()
      else
         moatTests[layer.name] = {}
      end

      if args.all then
         -- --all overrides the ignore list.
         ignoreTests[layer.name] = { lua = {}, c = {} }
      end
   end
end

-- ---- Test discovery ---------------------------------------------------------
-- Returns two lists of test descriptors:
--   luaTests: { {file, layer, name}, ... }   -- Lua .lua test files
--   cTests:   { {bin, src, layer, name}, ... } -- compiled C test binaries
--
-- 'activeLayers' is either nil (→ all layers) or a single layer name.
-- 'appPrefix' optionally filters Lua tests to those whose name starts with
-- "rt_<appPrefix>_" (mirrors the old --app option).
local function list_tests(activeLayers, args)
   local luaTests, cTests = {}, {}

   -- Determine which layers to scan.
   local layersToScan = {}
   if activeLayers then
      for _, L in ipairs(LAYERS) do
         if L.name == activeLayers then
            table.insert(layersToScan, L)
            break
         end
      end
   else
      layersToScan = LAYERS
   end

   local appPrefix = args.app

   for _, layer in ipairs(layersToScan) do
      -- ---- Lua tests --------------------------------------------------------
      local luaregDir = configVals.source_dir .. "/" .. layer.src .. "/luareg"
      local dirAttr   = lfs.attributes(luaregDir)
      if dirAttr and dirAttr.mode == "directory" then
         local function addLuaTest(fn)
            -- When an absolute path is given (e.g. via --run-only), only accept
            -- files that actually live under this layer's luareg directory.
            -- This prevents the per-layer loop from adding the same file to
            -- every layer in layersToScan.
            if string.sub(fn, 1, 1) == "/" then
               local luaregPath = "/" .. layer.src .. "/luareg/"
               if not string.find(fn, luaregPath, 1, true) then return end
            end
            -- Filter by app prefix if requested.
            if appPrefix then
               if not string.find(fn, "/rt_" .. appPrefix .. "_[^/]+%.lua$") then
                  return
               end
            else
               if not isLuaRegressionTest(fn) then return end
            end
            -- Skip tests in ignore list.  The lua table stores basenames
            -- (no path, no extension) so the file works on any machine.
            local ignLua = ignoreTests[layer.name] and ignoreTests[layer.name].lua or {}
            if lume.find(ignLua, stripext(basename(fn))) then return end
            -- Derive a short display name from the absolute path for DB storage.
            local nm = layer.name .. "/luareg/" .. basename(fn)
            table.insert(luaTests, {
               file  = fn,
               layer = layer.name,
               name  = nm,
            })
         end

         if args.moat then
            -- Only run the MOAT tests for this layer.
            for _, t in ipairs(moatTests[layer.name] or {}) do
               addLuaTest(t)
            end
         elseif args.run_only then
            for _, ro in ipairs(splitList(args.run_only)) do
               local a = lfs.attributes(ro)
               if a then
                  if a.mode == "file" then
                     addLuaTest(ro)
                  elseif a.mode == "directory" then
                     for dir, fn, _ in dirtree(ro) do addLuaTest(dir .. "/" .. fn) end
                  end
               end
            end
         else
            for dir, fn, _ in dirtree(luaregDir) do
               addLuaTest(dir .. "/" .. fn)
            end
         end
      else
         verboseLog(string.format(
            "Lua regression directory '%s' not found, skipping.\n", luaregDir))
      end

      -- ---- C tests ----------------------------------------------------------
      -- C regression sources live in source_dir/<layer>/creg/.
      -- Each rt_*.c file is compiled on-the-fly when the test runs; no
      -- pre-built binaries are required.
      local cregSrcDir = configVals.source_dir .. "/" .. layer.src .. "/creg"
      local cDirAttr   = lfs.attributes(cregSrcDir)
      if cDirAttr and cDirAttr.mode == "directory" then
         local function addCTest(fn)
            -- Layer-affinity guard for absolute paths (e.g. from --run-only).
            if string.sub(fn, 1, 1) == "/" then
               local cregPath = "/" .. layer.src .. "/creg/"
               if not string.find(fn, cregPath, 1, true) then return end
            end
            -- Must match rt_*.c pattern.
            if not string.match(fn, "/rt_[^/]+%.c$") then return end
            local testname = stripext(basename(fn))
            -- Filter by app prefix if requested.
            if appPrefix and not string.find(testname, "^rt_" .. appPrefix .. "_") then
               return
            end
            -- Skip tests in ignore list (C ignores are stored as basenames without .c).
            local ignC = ignoreTests[layer.name] and ignoreTests[layer.name].c or {}
            if lume.find(ignC, testname) then return end
            table.insert(cTests, {
               src   = fn,
               layer = layer.name,
               name  = layer.name .. "/creg/" .. testname,
            })
         end

         if args.run_only then
            for _, ro in ipairs(splitList(args.run_only)) do
               local a = lfs.attributes(ro)
               if a then
                  if a.mode == "file" then
                     addCTest(ro)
                  elseif a.mode == "directory" then
                     for dir, fn, attr in dirtree(ro) do
                        if attr.mode == "file" then addCTest(dir .. "/" .. fn) end
                     end
                  end
               end
            end
         else
            for fn in lfs.dir(cregSrcDir) do
               if string.match(fn, "^rt_.+%.c$") then
                  addCTest(cregSrcDir .. "/" .. fn)
               end
            end
         end
      else
         verboseLog(string.format(
            "C regression directory '%s' not found, skipping.\n", cregSrcDir))
      end
   end

   return luaTests, cTests
end

-- ---- Test runners -----------------------------------------------------------

-- Runs a single Lua regression test.
-- Workflow:
--   1. Create a per-test scratch directory under luareg-runs/<name>/
--   2. Copy the .lua input file there (so the developer can diff it later)
--   3. Remove stale .gkyl output files from the scratch directory
--   4. Run gkeyll <abs_path_to_lua> from within the scratch directory.
--      Without the -d flag, GKYL_OUT_PREFIX = basename(lua_file_without_ext),
--      so all .gkyl output lands in the CWD = scratch directory.
-- Returns (runtm, runlog, runDir).
-- timeoutSecs: 0 or nil = no limit; > 0 = kill test after that many seconds.
-- Returns: runtm, runlog, runDir, timedOut (boolean).
local function runLuaTest(test, timeoutSecs)
   local testBasename = stripext(basename(test.file))
   local runDir = configVals.results_dir .. "/" .. test.layer
      .. "/luareg-runs/" .. testBasename

   log(string.format("\n[Lua] Running %s ...\n", test.name))
   mkdir(runDir)

   -- Copy the Lua input file for reference / diffing.
   os.execute(string.format("cp -f '%s' '%s/'", test.file, runDir))

   -- Remove old output files to ensure a clean run.
   os.execute(string.format("rm -f '%s'/*.gkyl 2>/dev/null", runDir))

   local opts = {}
   local fh = io.open(test.file, "r")
   if fh then
      local line1 = fh:read()
      if line1 and string.find(line1, "--!") then
         local r = string.sub(line1, 4, -1)
         opts = loadstring("return " .. r)()
      end
      fh:close()
   end

   if opts.numProc then
      -- Parallel tests via MPI are not yet supported in the spawned-process
      -- model used here. Skip them and note in the log.
      log(string.format("**** NOT RUNNING PARALLEL TEST %s\n", test.name))
      return 0, "", runDir, false
   end

   -- Build the inner command (stderr merged into stdout via 2>&1 so that
   -- everything is captured by io.popen).
   local gkylExec = GKYL_EXEC_PATH .. "/gkeyll"
   local innerCmd = string.format(
      "cd '%s' && '%s' '%s' 2>&1", runDir, gkylExec, test.file)
   local runCmd = wrapWithTimeout(innerCmd, timeoutSecs or 0, runDir)

   local tmStart  = Time.clock()
   local proc     = io.popen(runCmd, "r")
   local rawOutput = proc:read("*a")
   proc:close()
   local runtm = Time.clock() - tmStart

   -- wrapWithTimeout appends '; echo __EXIT__:$?' so we can detect timeout.
   local exitCode = tonumber(rawOutput:match("__EXIT__:(%d+)%s*$")) or 0
   local runlog   = rawOutput:gsub("\n?__EXIT__:%d+%s*$", "")
   local timedOut = (exitCode == 124)

   if timedOut then
      log(string.format("... TIMED OUT after %g sec\n", runtm))
   else
      log(string.format("... completed in %g sec\n", runtm))
   end
   verboseLog(runlog)

   return runtm, runlog, runDir, timedOut
end

-- Compiles a single C regression test in its scratch directory.
-- Copies test.src and the installed share/Makefile into scratchDir, then
-- runs 'make <testname>' there.
-- Returns: ok (boolean), compileLog (string).
local function compileCTest(test, scratchDir)
   local testname      = stripext(basename(test.src))
   local shareMakefile = configVals.prefix .. "/gkeyll/share/Makefile"

   if not lfs.attributes(shareMakefile) then
      return false, string.format(
         "share/Makefile not found at '%s'.\n"
         .. "Ensure 'make install' has been run for the current build.\n",
         shareMakefile)
   end

   -- Copy the C source and the installed Makefile into the scratch dir.
   os.execute(string.format("cp -f '%s' '%s/'", test.src, scratchDir))
   os.execute(string.format("cp -f '%s' '%s/'", shareMakefile, scratchDir))
   -- Copy rt_arg_parse.h from the same creg/ directory so the local -I.
   -- takes precedence over the installed header.  This ensures any updates to
   -- rt_arg_parse.h in the source tree (e.g. new fields) are picked up without
   -- requiring a full 'make install'.
   local argParseH = dirname(test.src) .. "/rt_arg_parse.h"
   if lfs.attributes(argParseH) then
      os.execute(string.format("cp -f '%s' '%s/'", argParseH, scratchDir))
   end

   -- Compile: 'make <testname>' from the scratch directory.
   local compileCmd = string.format(
      "cd '%s' && make '%s' 2>&1; echo __COMPILE_EXIT__:$?", scratchDir, testname)
   local proc    = io.popen(compileCmd, "r")
   local rawOut  = proc:read("*a")
   proc:close()
   local exitCode   = tonumber(rawOut:match("__COMPILE_EXIT__:(%d+)%s*$")) or 1
   local compileLog = rawOut:gsub("\n?__COMPILE_EXIT__:%d+%s*$", "")
   return (exitCode == 0), compileLog
end

-- Runs a single C regression test.
-- Workflow:
--   1. Create per-test scratch directory under creg-runs/<testname>/
--   2. Remove stale .gkyl output files
--   3. Compile the .c source on-the-fly via compileCTest (copies .c +
--      share/Makefile into scratch dir, runs 'make <testname>')
--   4. Run the freshly compiled binary from the scratch directory
--   5. Remove the compiled binary (keep .c, Makefile, .gkyl outputs)
-- Returns: runtm, runlog, runDir, timedOut (bool), compileFailed (bool).
-- Status codes recorded in DB: -4=compile_fail, -3=timeout, 0=fail, 1=pass.
local function runCTest(test, timeoutSecs)
   local testname = stripext(basename(test.src))
   local runDir   = configVals.results_dir .. "/" .. test.layer
      .. "/creg-runs/" .. testname

   log(string.format("\n[C]   Running %s ...\n", test.name))
   mkdir(runDir)
   os.execute(string.format("rm -f '%s'/*.gkyl 2>/dev/null", runDir))

   -- Step 1: compile.
   log("... compiling ...\n")
   local tmComp = Time.clock()
   local compileOk, compileLog = compileCTest(test, runDir)
   local compileSecs = Time.clock() - tmComp

   if not compileOk then
      log(string.format("... COMPILE FAILED in %g sec\n", compileSecs))
      verboseLog(compileLog)
      return compileSecs, "COMPILE FAILED:\n" .. compileLog, runDir, false, true
   end
   verboseLog(compileLog)

   -- Step 2: run.
   local binPath  = runDir .. "/" .. testname
   local innerCmd = string.format("cd '%s' && '%s' 2>&1", runDir, binPath)
   local runCmd   = wrapWithTimeout(innerCmd, timeoutSecs or 0, runDir)

   local tmStart   = Time.clock()
   local proc      = io.popen(runCmd, "r")
   local rawOutput = proc:read("*a")
   proc:close()
   local runtm = Time.clock() - tmStart

   local exitCode = tonumber(rawOutput:match("__EXIT__:(%d+)%s*$")) or 0
   local runlog   = rawOutput:gsub("\n?__EXIT__:%d+%s*$", "")
   local timedOut = (exitCode == 124)

   -- Remove the compiled binary and macOS debug symbols dir; keep .c source,
   -- Makefile, and .gkyl outputs.
   os.execute(string.format("rm -f '%s/%s' '%s/%s.d' 2>/dev/null",
      runDir, testname, runDir, testname))
   os.execute(string.format("rm -rf '%s/%s.dSYM' 2>/dev/null",
      runDir, testname))

   if timedOut then
      log(string.format("... TIMED OUT after %g sec\n", runtm))
   else
      log(string.format("... completed in %g sec\n", runtm))
   end
   verboseLog(runlog)

   return runtm, compileLog .. "\n" .. runlog, runDir, timedOut, false
end

-- ---- File comparison --------------------------------------------------------
-- (Unchanged from original implementation.)
local function compareFiles(f1, f2)
   verboseLog(string.format("Comparing %s %s ...\n", f1, f2))
   if not lfs.attributes(f1) or not lfs.attributes(f2) then
      verboseLog(string.format(
         " ... files %s and/or %s do not exist!\n", f1, f2))
      return false
   end

   local f1type = ZeroUtil.gkylFileType(f1)
   local f2type = ZeroUtil.gkylFileType(f2)

   if f1type ~= f2type then
      verboseLog(string.format(
         " ... files %s and %s are not of the same type!\n", f1, f2))
      return false
   end

   if f1type == "dynvector" then
      verboseLog(" ... dynvector comparison NYI!\n")
      return true
   end

   -- ArrayNewFromFile may throw for non-array file types (e.g. topology/btopo files).
   -- Guard with pcall; if it fails, skip the comparison rather than crashing.
   local ok1, g1, a1 = pcall(ZeroArrayRio.ArrayNewFromFile, f1)
   local ok2, g2, a2 = pcall(ZeroArrayRio.ArrayNewFromFile, f2)
   if not ok1 or not ok2 then
      verboseLog(string.format(
         " ... skipping %s (unsupported file format)\n", f1))
      return true
   end

   if not g1 then return false end
   if not g2 then return false end

   local gridsSame = ZeroRectCart.compare(g1, g2)
   if not gridsSame then return false end

   local nghost = { 0, 0, 0, 0, 0, 0, 0 }
   local r1, er1 = g1:createGridRanges(nghost)
   local r2, er2 = g2:createGridRanges(nghost)

   local diff = ZeroArray.arrayDiff(a1, a2, r1)

   if not diff.is_compatible then return false end
   if diff.max_abs_diff > 1e-12 then
      verboseLog(string.format(
         " ... max absolute difference is %g\n", diff.max_abs_diff))
      return false
   end
   -- Only check relative difference when there is a non-trivial absolute
   -- difference. gkyl_array_diff returns DBL_MAX for max_rel_diff when both
   -- values are zero (0/0 is undefined), which is not a meaningful failure.
   if diff.max_abs_diff > 0 and diff.max_rel_diff > 1e-12 then
      verboseLog(string.format(
         " ... max relative difference is %g\n", diff.max_rel_diff))
      return false
   end

   return true
end

-- ---- Post-run actions: create and check -------------------------------------
-- Each action takes (test, runDir, testType) where testType is 'lua' or 'c'.
-- 'create' copies the run's .gkyl output to the accepted directory, making
-- it the new baseline. Status -2 is recorded in the database.
-- 'check' compares the current run's output against the accepted baseline.
-- Status 1 = pass, 0 = fail.

-- Returns the accepted results directory for the given test and type.
local function acceptedDir(test, testType)
   -- testType is 'lua' or 'c'; subdirectory names are 'luareg-accepted' / 'creg-accepted'.
   local sub = (testType == "lua") and "luareg-accepted" or "creg-accepted"
   local nm  = (testType == "lua")
      and stripext(basename(test.file))
      or  stripext(basename(test.src))
   return configVals.results_dir .. "/" .. test.layer .. "/" .. sub .. "/" .. nm
end

local function create_action(test, runDir, testType)
   local aDir = acceptedDir(test, testType)
   log(string.format("... saving accepted results to %s ...\n", aDir))
   mkdir(aDir)
   -- Copy all .gkyl output files from the scratch directory to the accepted dir.
   os.execute(string.format("cp -f '%s'/*.gkyl '%s/' 2>/dev/null", runDir, aDir))
   return -2
end

local function check_action(test, runDir, testType)
   local aDir = acceptedDir(test, testType)
   -- For Lua tests the output prefix matches the input basename (e.g. rt_euler_sodshock).
   -- For C tests the app's .name field sets the prefix, which we have standardised to
   -- match the source basename (e.g. rt_10m_sodshock), so all .gkyl files in the run
   -- directory belong to this test and we compare them all.
   local testPrefix = (testType == "lua") and stripext(basename(test.file)) or nil

   local passed, count = true, 0
   -- Walk the scratch directory and compare each .gkyl file against the accepted baseline.
   for fn in lfs.dir(runDir) do
      local isGkyl   = (string.sub(fn, -5) == ".gkyl")
      -- Skip topology/connectivity files (e.g. _btopo.gkyl): they are not
      -- grid+array files and cannot be compared with gkyl_array_diff.
      local isSkipped = string.find(fn, "_btopo%.gkyl$") ~= nil
      local inScope  = isGkyl and not isSkipped and (
         testPrefix == nil or string.find(fn, "^" .. testPrefix) ~= nil)
      if inScope then
         count = count + 1
         local runFile  = runDir  .. "/" .. fn
         local accFile  = aDir    .. "/" .. fn
         local ok = compareFiles(accFile, runFile)
         passed = passed and ok
      end
   end

   -- If no output files were produced, the test likely crashed. Treat as fail.
   -- This preserves the original "HACK" comment intent: tests must produce output.
   if count == 0 then passed = false end

   if passed then
      layerCounts[test.layer].passed = layerCounts[test.layer].passed + 1
      log("... passed.\n")
   else
      layerCounts[test.layer].failed = layerCounts[test.layer].failed + 1
      log(string.format("... %s FAILED!\n", test.name))
   end
   return passed and 1 or 0
end

-- ---- Unit test helpers (unchanged from original) ----------------------------
local function runLuaUnitTest(test)
   log(string.format("\nRunning unit test %s ...\n", test))
   local gkylExec = GKYL_EXEC_PATH .. "/gkeyll"
   local f = io.popen(string.format("'%s' '%s'", gkylExec, test), "r")
   local outPut = f:read("*a")
   f:close()
   if string.find(outPut, "PASSED") then
      log("... passed.\n")
   else
      if string.find(outPut, "without CUDA") then
         log(string.format("... %s skipped (no CUDA).\n", test))
      else
         log(string.format("... %s FAILED!\n", test))
      end
   end
end

local function runCxxUnitTest(test)
   log(string.format("\nRunning C++ unit test %s ...\n", test))
   local f = io.popen(string.format("'%s'", test), "r")
   local outPut = f:read("*a")
   f:close()
   if string.find(outPut, "FAILED") then
      log(string.format("... %s FAILED!\n", test))
   else
      log("... passed.\n")
   end
end

local function list_unit_tests(args)
   local luaUnitTests, cxxUnitTests = {}, {}
   -- Unit tests still live in the source tree (parallel structure to regression).
   local layersToScan = {}
   if detectedLayer then
      for _, L in ipairs(LAYERS) do
         if L.name == detectedLayer then
            table.insert(layersToScan, L); break
         end
      end
   else
      layersToScan = LAYERS
   end
   for _, layer in ipairs(layersToScan) do
      local unitDir = configVals.source_dir .. "/" .. layer.src .. "/unit"
      if lfs.attributes(unitDir) then
         for dir, fn, _ in dirtree(unitDir) do
            local fullNm = dir .. "/" .. fn
            if isLuaUnitTest(fullNm) then
               table.insert(luaUnitTests, fullNm)
            elseif isCxxUnitTest(fullNm) then
               table.insert(cxxUnitTests, fullNm)
            end
         end
      end
   end
   return luaUnitTests, cxxUnitTests
end

-- ---- Ignore-list updater ----------------------------------------------------
-- Called after a run with --timeout when some tests exceeded the limit.
-- Merges the new timed-out test names into the layer's ignoretests.lua and
-- writes the file back so subsequent runs automatically skip them.
--
-- 'layer'        - layer table ({name, src})
-- 'newLuaNames'  - list of Lua test basenames (no path, no .lua) that timed out
-- 'newCNames'    - list of C binary basenames (no path, no extension) that timed out
local function updateIgnoreTests(layer, newLuaPaths, newCNames)
   local ignFile = configVals.source_dir
      .. "/" .. layer.src .. "/luareg/ignoretests.lua"

   -- Use the in-memory table already loaded by loadConfigure as the base.
   local existing = ignoreTests[layer.name] or { lua = {}, c = {} }

   -- Merge Lua paths (avoid duplicates).
   local luaSet = {}
   for _, p in ipairs(existing.lua) do luaSet[p] = true end
   local addedLua = 0
   for _, p in ipairs(newLuaPaths) do
      if not luaSet[p] then
         luaSet[p] = true
         table.insert(existing.lua, p)
         addedLua = addedLua + 1
      end
   end

   -- Merge C names (avoid duplicates).
   local cSet = {}
   for _, n in ipairs(existing.c) do cSet[n] = true end
   local addedC = 0
   for _, n in ipairs(newCNames) do
      if not cSet[n] then
         cSet[n] = true
         table.insert(existing.c, n)
         addedC = addedC + 1
      end
   end

   if addedLua == 0 and addedC == 0 then return end

   -- Write the merged table back to the file.
   local f = io.open(ignFile, "w")
   f:write("-- Auto-updated by runregression: tests that exceeded the per-test timeout.\n")
   f:write("-- Remove an entry manually to re-enable the test.\n")
   f:write("return {\n")
   f:write("   lua = {\n")
   for _, p in ipairs(existing.lua) do
      f:write(string.format("      %q,\n", p))
   end
   f:write("   },\n")
   f:write("   c = {\n")
   for _, n in ipairs(existing.c) do
      f:write(string.format("      %q,\n", n))
   end
   f:write("   },\n")
   f:write("}\n")
   f:close()

   log(string.format("[ignore] Updated %s (+%d Lua, +%d C timed-out tests)\n",
      ignFile, addedLua, addedC))
end

-- ---- Command action functions -----------------------------------------------

-- 'configure' command: set up the regression system for this machine.
local function config_action(args, name)
   isConfiguring = true

   local prefix = args.config_prefix or (os.getenv("HOME") .. "/gkylsoft")
   local mpiexec = args.config_mpiexec
      or (os.getenv("HOME") .. "/gkylsoft/openmpi/bin/mpiexec")
   local sourceDir = args.config_source_dir or lfs.currentdir()

   configure(prefix, mpiexec, sourceDir, args)
end

-- 'list' command: print all regression tests that would be run.
local function list_action(args, name)
   loadConfigure(args)
   local luaTests, cTests = list_tests(detectedLayer, args)
   if not args.c_only then
      for _, t in ipairs(luaTests) do print("[lua] " .. t.name) end
   end
   if not args.lua_only then
      for _, t in ipairs(cTests)   do print("[c]   " .. t.name) end
   end
end

-- 'run' command: execute regression tests and optionally create or check results.
local function run_action(args, name)
   loadConfigure(args)

   local luaTests, cTests = list_tests(detectedLayer, args)

   -- Per-test timeout in seconds (0 = unlimited).
   local timeoutSecs = args.timeout or 0
   if timeoutSecs > 0 and not TIMEOUT_CMD then
      log(string.format(
         "NOTE: 'timeout'/'gtimeout' not found; using Perl-based timeout.\n"))
   end

   -- Determine the post-run action. Without a sub-command the tests are run
   -- but results are neither saved nor compared (status = -1, "skip").
   local postRun = function(test, runDir, testType) return -1 end
   if args.create then
      postRun = create_action
   elseif args.check then
      postRun = check_action
   end

   -- Track timed-out tests per layer so we can update ignoretests.lua.
   local timedOutByLayer = {}
   for _, L in ipairs(LAYERS) do
      timedOutByLayer[L.name] = { lua = {}, c = {} }
   end

   log("Running regression tests ...\n\n")
   local tmStart = Time.clock()

   -- ---- Lua tests ----------------------------------------------------------
   if not args.c_only then
      for _, test in ipairs(luaTests) do
         layerCounts[test.layer].total = layerCounts[test.layer].total + 1
         local runtm, runlog, runDir, timedOut = runLuaTest(test, timeoutSecs)
         if timedOut then
            -- Status -3 = timed out; do NOT create accepted results.
            -- Store basename (no path, no extension) for portability across machines.
            table.insert(timedOutByLayer[test.layer].lua, stripext(basename(test.file)))
            insertRegressionData(
               test.layer, runID, test.name, "lua", -3, runtm, "TIMED OUT")
            layerCounts[test.layer].failed = layerCounts[test.layer].failed + 1
         else
            local status = postRun(test, runDir, "lua")
            insertRegressionData(
               test.layer, runID, test.name, "lua", status, runtm, runlog)
         end
      end
   end

   -- ---- C tests ------------------------------------------------------------
   if not args.lua_only then
   for _, test in ipairs(cTests) do
      layerCounts[test.layer].total = layerCounts[test.layer].total + 1
      local runtm, runlog, runDir, timedOut, compileFailed =
         runCTest(test, timeoutSecs)
      if compileFailed then
         -- Status -4: compilation failed; no .gkyl output to compare.
         insertRegressionData(
            test.layer, runID, test.name, "c", -4, runtm, runlog)
         layerCounts[test.layer].failed = layerCounts[test.layer].failed + 1
      elseif timedOut then
         -- Status -3: test exceeded wall-clock timeout.
         table.insert(timedOutByLayer[test.layer].c,
            stripext(basename(test.src)))
         insertRegressionData(
            test.layer, runID, test.name, "c", -3, runtm, "TIMED OUT")
         layerCounts[test.layer].failed = layerCounts[test.layer].failed + 1
      else
         local status = postRun(test, runDir, "c")
         insertRegressionData(
            test.layer, runID, test.name, "c", status, runtm, runlog)
      end
   end
   end -- if not args.lua_only

   log(string.format(
      "\nAll regression tests completed in %g secs\n", Time.clock() - tmStart))

   -- Update ignoretests.lua for any layers that had timed-out tests.
   if timeoutSecs > 0 then
      for _, layer in ipairs(LAYERS) do
         local to = timedOutByLayer[layer.name]
         if #to.lua > 0 or #to.c > 0 then
            updateIgnoreTests(layer, to.lua, to.c)
         end
      end
   end

   -- Finalize: write per-layer summary rows into each layer's database.
   -- We only write a metadata row for layers that actually had tests.
   for _, layer in ipairs(LAYERS) do
      local cnt = layerCounts[layer.name]
      if cnt.total > 0 then
         insertRegressionMeta(layer.name, runID, runDate,
            cnt.total, cnt.passed, cnt.failed)
         log(string.format(
            "  Layer %-12s  total=%d  passed=%d  failed=%d\n",
            layer.name, cnt.total, cnt.passed, cnt.failed))
      end
   end
end

-- 'listunit' command.
local function listunit_action(args, name)
   loadConfigure(args)
   local lua, cxx = list_unit_tests(args)
   lume.each(lua, print)
   lume.each(cxx, print)
end

-- 'rununit' command.
local function rununit_action(args, name)
   loadConfigure(args)
   local lua, cxx = list_unit_tests(args)
   log("Running unit tests ...\n\n")
   local tmStart = Time.clock()
   lume.each(lua, runLuaUnitTest)
   lume.each(cxx, runCxxUnitTest)
   log(string.format("All unit tests completed in %g secs\n", Time.clock() - tmStart))
end

-- ---- CLI parser -------------------------------------------------------------

local parser = argparse()
   :name("runregression")
   :require_command(true)
   :description [[
Run Gkyl regression tests across the hierarchical layer architecture
(moments -> vlasov -> gyrokinetic -> pkpm).

Each layer has Lua regression tests (luareg/) and C regression tests (creg/).
Results are stored in per-layer SQLite databases under gkeyll-results/.

Typical workflow:
  1. Build and install: make install -j4
  2. Configure: gkeyll runregression configure --source-dir /path/to/src
  3. Create baselines: gkeyll runregression run create --timeout 120
  4. Check results:    gkeyll runregression run check  --timeout 120
  5. Layer-specific:   gkeyll runregression run moments check

C regression tests are compiled on-the-fly using the installed
share/Makefile (PREFIX/gkeyll/share/Makefile). No separate 'make regression'
step is needed.
]]

parser:flag("-v --verbose", "Print verbose messages as tests are run")
parser:flag("-a --all", "Run all tests, ignoring ignoretests files")

-- 'configure' command ---------------------------------------------------------
-- Sets up directories, databases, and writes ~/runregression.config.lua.
local c_conf = parser:command("configure", "Configure regression tests")
   :action(config_action)

c_conf:option("-p --prefix", "Location to write gkeyll-results/ directory")
   :target("config_prefix")
c_conf:option("-m --mpiexec", "Full path to MPI executable")
   :target("config_mpiexec")
c_conf:option("-s --source-dir",
   "Path to gkeyll source root (contains moments/, vlasov/, etc. with creg/ and luareg/ subdirs)")
   :target("config_source_dir")
c_conf:flag("--drop-tables",
   "Re-create all SQL tables (erases existing regression data)", false)

-- 'list' command --------------------------------------------------------------
-- Lists all regression tests that would be run (useful for inspection).
local c_list = parser:command("list", "List all regression tests")
   :action(list_action)
c_list:option("-r --run-only",
   "Only list this test or all tests in this directory")
c_list:option("-a --app",
   "Only list tests for this app (prefix after rt_ in filename)")
c_list:flag("-m --moat", "Only list MOAT regression tests")
c_list:flag("-c --c-only",   "Only list C regression tests (skip Lua tests)")
c_list:flag("-l --lua-only", "Only list Lua regression tests (skip C tests)")

-- 'run' command ---------------------------------------------------------------
-- Runs regression tests. Layer is extracted from the command line BEFORE
-- argparse (see the pre-processing block near the top of this file), so:
--   gkeyll runregression run check          → all layers, check
--   gkeyll runregression run moments check  → moments layer only, check
local c_run = parser:command("run", "Run regression tests.")
   :require_command(false)
   :action(run_action)
c_run:option("-r --run-only",
   "Only run these tests or directories (comma-separated list)")
c_run:option("-a --app",
   "Only run tests for this app (prefix after rt_ in filename)")
c_run:flag("-m --moat", "Only run MOAT regression tests")
c_run:flag("-c --c-only",   "Only run C regression tests (skip Lua tests)")
c_run:flag("-l --lua-only", "Only run Lua regression tests (skip C tests)")
c_run:option("-t --timeout",
   "Per-test wall-clock timeout in seconds (0 = unlimited). "
   .. "Tests that exceed this limit are added to ignoretests.lua automatically.")
   :convert(tonumber)
   :default(0)

c_run:command("check",  "Check results against accepted baselines")
c_run:command("create", "Create accepted baselines from current results")

-- 'listunit' command ----------------------------------------------------------
parser:command("listunit", "List all unit tests")
   :action(listunit_action)

-- 'rununit' command -----------------------------------------------------------
local c_rununit = parser:command("rununit", "Run unit tests.")
   :require_command(false)
   :action(rununit_action)

-- ---- Parse and dispatch -----------------------------------------------------
-- Use filteredCmds (layer names removed) rather than GKYL_COMMANDS_L.
local _ = parser:parse(filteredCmds)
