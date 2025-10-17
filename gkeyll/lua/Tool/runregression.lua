-- Gkyl ------------------------------------------------------------------------
--
-- Runs regression tests and compares with accepted results.
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

local currDir = lfs.currentdir()
if not string.match(currDir, "regression$") then
   print("Can only run from regression source directory.")
   return 1
end

local Alloc = require "Lib.Alloc"
local Logger = require "Lib.Logger"
local Time = require "Lib.Time"
local argparse = require "Lib.argparse"
local date = require "xsys.date"
local lume = require "Lib.lume"
local Lin = require "Lib.Linalg"
local sql = require "sqlite3"

local ZeroArray = require "DataStruct.ZeroArray"
local ZeroArrayRio = require "DataStruct.ZeroArrayRio"
local ZeroDynVector = require "DataStruct.ZeroDynVector"
local ZeroRectCart = require "Grid.ZeroRectCart"
local ZeroUtil = require "Lib.ZeroUtil"

-- need to change this as it is keyed on input file name
GKYL_OUT_PREFIX = lfs.currentdir() .. "/" .. "runregression"

local log = Logger { logToFile = true }
local verboseLog = function (msg) end -- default no messages are written
local verboseLogger = function (msg) log(msg) end

-- number of passed/failed and total regression tests
local numFailedTests = 0
local numPassedTests = 0
local numTotalTests = 0

-- number of passed/failed and total unit tests
local numFailedUnitTests = 0
local numPassedUnitTests = 0
local numTotalUnitTests = 0

-- flag to indicate if we are configuring system
local isConfiguring = false
-- global value to store config information
local configVals = nil
-- global list of tests to ignore
local ignoreTests = {}
-- list of MOAT regression
local moatTests = {}
-- Sqlite3 database for storing regression data
local sqlConn = nil
-- stored procedures to store data in tables
local insertRegressionDataProc = nil
local insertRegressionMetaProc = nil

-- Name of configuration file
local confFile = os.getenv("HOME") .. "/runregression.config.lua"

-- function to split comma separated list (string) into table
local function splitList(listStr)
   local words = {}
   for w in listStr:gmatch('[^,%s]+') do
      table.insert(words, w)
   end
   return words
end

-- generate ID for this run of regression system
local runID = uuid()
-- date when tests were run
local runDate = date(false):fmt("${iso}")

-- Table structure for data stored in SQLite DB
-- 
-- table RegressionMeta (
--   guid text,
--   tstamp text,
--   GKYL_EXEC text,
--   GKYL_HG_CHANGESET text,
--   GKYL_BUILD_DATE text,
--   ntotal integer,
--   npass integer,
--   nfail integer
-- );
--
-- table RegressionData (
--   guid text,
--   name text,
--   status integer,
--   runtime real,
--   runlog text
-- );
--

-- loads configuration file
local function loadConfigure(args)
   local f = loadfile(confFile)
   if not f then
      log("Regression tests not configured! Run config command first\n")
      os.exit(1)
   end
   configVals = f()

   if args.verbose then
      verboseLog = verboseLogger -- set verbose logger
   end

   local g = loadfile("ignoretests.lua")
   if not args.all then
      ignoreTests = g()
   end

   local m = loadfile("moat.lua")
   moatTests = m()

   -- open connection to SQL DB
   sqlConn = sql.open(string.format("%s/regressiondb", configVals['results_dir']))
   -- SQL stored procedures to insert data into table
   insertRegressionDataProc = sqlConn:prepare [[
     insert into RegressionData values (
       ?, ?, ?, ?, ?
     )
   ]]
   insertRegressionMetaProc = sqlConn:prepare [[
     insert into RegressionMeta values (
       ?, ?, ?, ?, ?, ?, ?, ?
     )
   ]]
   
   return configVals
end

-- functions to insert data into the tables
local function insertRegressionData(id, nm, status, runtm, runlog)
   insertRegressionDataProc:reset():bind(
      id,
      nm,
      status,
      runtm,
      runlog):step()
end
local function insertRegressionMeta(id, tm, ntotal, npass, nfail)
   insertRegressionMetaProc:reset():bind(
      id,
      tm,
      GKYL_EXEC,
      GKYL_GIT_CHANGESET,
      GKYL_BUILD_DATE,      
      ntotal,
      npass,
      nfail):step()
end

-- creates directories if they don't exist
local function mkdir(path)
   local sep, pStr = package.config:sub(1,1), "/"
   for dir in path:gmatch("[^" .. sep .. "]+") do
      pStr = pStr .. dir .. sep
      lfs.mkdir(pStr)
   end
end

-- configure paths
local function configure(prefix, mpiExec, args)
   local mpiAttr = lfs.attributes(mpiExec)

   -- FOR NOW THIS IS DISABLED TILL I FIGURE OUT HOW TO DO PARALLEL
   -- TESTS (AH, 5/24/2019)
   if false then
      if mpiAttr and mpiAttr.mode == "file" and string.sub(mpiAttr.permissions, 3, 3) == "x" then
	 log(string.format("Setting MPIEXEC to %s ...\n", mpiExec))
      else
	 assert(false, string.format("MPI launch binary %s not found or is not an executable!", mpiExec))
      end
   end

   local prefixAttr = lfs.attributes(prefix)
   if prefixAttr and prefixAttr.mode == "directory" then
      log(string.format("Will write accepted results to %s/gkyl-results ...\n", prefix))
   else
      assert(false, string.format("Prefix %s is not a directory!", prefix))
   end

   -- create directory (only if it does not exist)
   mkdir(string.format("%s/gkyl-results", prefix))

   -- regression data is stored in SQLite. Create tables if needed
   local regressionDb = string.format("%s/gkyl-results/regressiondb", prefix)
   local regressionDbAttr = lfs.attributes(regressionDb)

   if args.drop_tables or regressionDbAttr==nil then
      log(string.format("Creating DB file %s\n", regressionDb))
      -- create Sqlite3 database to store regression data
      local conn = sql.open(regressionDb)
      
      -- Following tables are created: RegressionMeta that stores date
      -- when test was run and how many passed or
      -- failed. RegressionData stores output from each test. GUIDs
      -- are shared between tables so one can get all tests and their
      -- data given the GUID
      conn:exec [[

      drop table if exists RegressionMeta;
      create table RegressionMeta (
        guid text,
        tstamp text,
        GKYL_EXEC text,
        GKYL_HG_CHANGESET text,
        GKYL_BUILD_DATE text,
        ntotal integer,
        npass integer,
        nfail integer
      );

      drop table if exists RegressionData;
      create table RegressionData (
        guid text,
        name text,
        status integer,
        runtime real,
        runlog text
      );
     ]]
   end

   -- write information into config file
   local fn = io.open(confFile, "w")
   fn:write("return {\n")
   fn:write(string.format("mpiExec = \"%s\",\n", mpiExec))
   fn:write(string.format("results_dir = \"%s/gkyl-results\",\n", prefix))
   fn:write("}")
   fn:close()
end

-- Walks down a directory tree recursively: returns a coroutine that
-- yields directory, file-name and attribute object
local function dirtree(dir)
   assert(dir and dir ~= "", "directory parameter is missing or empty")
   if string.sub(dir, -1) == "/" then
      dir = string.sub(dir, 1, -2)
   end

   local function yieldtree(dir)
      for fn in lfs.dir(dir) do
	 if fn ~= "." and fn ~= ".." then
	    local fullNm = dir .. "/" .. fn
	    local attr = lfs.attributes(fullNm)
	    coroutine.yield(dir, fn, attr)
	    if attr.mode == "directory" then
	       yieldtree(fullNm)
	    end
	 end
      end
   end
   return coroutine.wrap(function() yieldtree(dir) end)
end

-- copies all files with specified extension to target directory
local function copyAllFiles(srcPath, ext, destPath)
   -- for now using os.execute to run the "cp" command. Perhaps this
   -- is not the best way to do this and one might want to use a more
   -- Lua-ish way
   local cp = "cp -rf" .. " " .. srcPath .. "-*." .. ext .. " " .. destPath
   verboseLog(string.format("... executing %s ...\n", cp))
   os.execute(cp)
end

-- runs a single lua test
local function runLuaTest(test)
   log(string.format("\nRunning test %s ...\n", test))

   numTotalTests = numTotalTests+1

   local opts = {}
   -- open Lua file and check if first line is a command line
   local f = io.open(test, "r")
   local line1 = f:read()
   if string.find(line1, "--!") then
      local r = string.sub(line1, 4, -1)
      opts = loadstring("return " .. r)()
   end
   f:close()

   -- before running tests delete all old output
   local rm = "rm -rf " .. string.sub(test, 1, -5) .. "-*.gkyl"
   os.execute(rm)

   -- construct command to run
   local runCmd = string.format("%s -d %s", GKYL_EXEC, test)
   if opts.numProc then
      runCmd = string.format("%s -n %d %s %s", configVals.mpiExec, opts.numProc, GKYL_EXEC, test)
   end

   -- run test
   local tmStart = Time.clock()

   local runlog = ""
   if opts.numProc then
      -- skip for now: NEED TO USE MPI_COMM_SPAWN FOR THIS, I THINK!
      log(string.format("**** NOT RUNNING PARALLEL TEST %s \n", runCmd))
   else
      local f = io.popen(runCmd, "r")
      for l in f:lines() do
	 if string.find(l, "*** LOAD ERROR ***") then
	    print("... input file error!")
	    print(f:read())
	 else
	    runlog = runlog .. l .. "\n"
	    verboseLog(l.."\n")
	 end
      end
   end
   local runtm = Time.clock()-tmStart
   log(string.format("... completed in %g sec\n", runtm))
   return runtm, runlog
end

-- runs a single shell-script test
local function runShellTest(test)
   log(string.format("Running shell test %s\n", test))
end

-- runs a single lua unit test
local function runLuaUnitTest(test)
   log(string.format("\nRunning test %s ...\n", test))

   local runCmd = string.format("%s %s", GKYL_EXEC, test)
   local f = io.popen(runCmd, "r")
   local outPut = f:read("*a")
   if string.find(outPut, "PASSED") then
      numPassedUnitTests = numPassedUnitTests+1
      log("... passed.\n")      
   else
      if string.find(outPut, "without CUDA") then
	 log(string.format(".... %s skipped as CUDA not present!\n", test))
      else
	 log(string.format("... %s FAILED!\n", test))
	 numFailedUnitTests = numFailedUnitTests+1
      end
   end
end

-- runs a single CXX unit test
local function runCxxUnitTest(test)
   log(string.format("\nRunning test %s ...\n", test))

   local runCmd = string.format("%s", test)
   local f = io.popen(runCmd, "r")
   local outPut = f:read("*a")
   if string.find(outPut, "FAILED") then
      log(string.format("... %s FAILED!\n", test))
      numFailedUnitTests = numFailedUnitTests+1
   else
      numPassedUnitTests = numPassedUnitTests+1
      log("... passed.\n")
   end
end

local function isLuaRegressionTest(fn, appPrefix)
   if appPrefix then
      return string.find(fn, "/rt_" .. appPrefix .. "_[^/]+%.lua$") ~= nil
   end
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

-- list of regression tests
local function list_tests(args)
   local luaRegTests, shellRegTests = {}, {}

   local appPrefix = args.app

   local function addTest(fn)
      if isLuaRegressionTest(fn, appPrefix) then
	 table.insert(luaRegTests, fn)
      end
   end

   if args.moat then
      -- only MOAT regressions are to be run
      for _, t in ipairs(moatTests) do addTest(t) end
   elseif args.run_only then
      -- get list of tests to run (this could be a comma separated
      -- list of file names or directory names, or combination of
      -- both)
      local runList = splitList(args.run_only)
      for _, ro in ipairs(runList) do
	 local a = lfs.attributes(ro)
	 if a.mode == "file" then
	    addTest(ro)
	 elseif a.mode == "directory" then
	    for dir, fn, _ in dirtree(ro) do addTest(dir .. "/" .. fn) end
	 end
      end
   else
      for dir, fn, _ in dirtree(".") do 
	 if (not dir:find("__needs_implementation")) then
	    addTest(dir .. "/" .. fn) 
	 end
      end
   end

   -- this function is used to filter out tests that should NOT be
   -- run. Note the confusing name
   local function shouldRun(t)
      if lume.find(ignoreTests, t) then return false end
      return true
   end
   
   return lume.filter(luaRegTests, shouldRun), shellRegTests
end

-- list of unit tests
local function list_unit_tests(args)
   local luaUnitTests, cxxUnitTests = {}, {}

   local function addTest(fn)
      if isLuaUnitTest(fn) then
	 table.insert(luaUnitTests, fn)
      elseif isCxxUnitTest(fn) then
	 table.insert(cxxUnitTests, fn)
      end
   end

   for dir, fn, _ in dirtree("../unit") do addTest(dir .. "/" .. fn) end

   return luaUnitTests, cxxUnitTests
end

-- function to handle "config" command
local function config_action(args, name)
   isConfiguring = true
   
   local prefix = args.config_prefix and args.config_prefix or
      os.getenv("HOME") .. "/gkylsoft"
   local mpiexec = args.config_mpiexec and  args.config_mpiexec or
      os.getenv("HOME") .. "/gkylsoft/openmpi/bin/mpiexec"

   configure(prefix, mpiexec, args)
end

-- function to handle "list" command
local function list_action(args, name)
   loadConfigure(args)
   
   local luaRegTests, shellRegTests = list_tests(args)
   lume.each(luaRegTests, print)
   lume.each(shellRegTests, print)
end

-- function to handle "create" sub-command of "run"
local function create_action(test)
   log(string.format("... creating accepted results for %s ...\n", test))
   local fullResultsDir = configVals.results_dir .. "/"
      .. string.sub(test, 1, -5) -- remove the initial ./ and last .lua

   verboseLog(string.format("... making directory %s ..\n", fullResultsDir))
   mkdir(fullResultsDir) -- recursively create all dirs as needed
   local srcPath = string.sub(test, 1, -5)
   copyAllFiles(srcPath, "gkyl", fullResultsDir)

   return -2
end

-- function to compare files
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
      verboseLog(string.format(" ... dynvector comparison NYI!"))
      return true
   end

   -- read arrays
   local g1, a1 = ZeroArrayRio.ArrayNewFromFile(f1)
   local g2, a2 = ZeroArrayRio.ArrayNewFromFile(f2)

   if not g1 then return false end
   if not g2 then return false end

   local gridsSame = ZeroRectCart.compare(g1, g2)
   if not gridsSame then
      return false
   end

   -- create range and extended ranges
   local nghost = { 0, 0, 0, 0, 0, 0, 0 }
   local r1, er1 = g1:createGridRanges(nghost)
   local r2, er2 = g2:createGridRanges(nghost)

   local diff = ZeroArray.arrayDiff(a1, a2, r1)

   -- check array differences
   if not diff.is_compatible then return false end
   if diff.max_abs_diff > 1e-12 then
      verboseLog(string.format(
		    " ... max absolute difference is %g\n", diff.max_abs_diff))
      return false
   end
   if diff.max_rel_diff > 1e-12 then
      verboseLog(string.format(
		    " ... max relative difference is %g\n", diff.max_rel_diff))
      return false
   end

   return true
end

-- function to handle "check" sub-command of "run"
local function check_action(test)
   local fullResultsDir = configVals.results_dir .. "/"
      .. string.sub(test, 1, -5) -- remove the initial ./ and last .lua

   local vloc = string.find(test, "/rt_[^/]+%.lua$")
   local outDirName = string.sub(test, 1, vloc)
   -- the very strange looking pattern below puts an escape character
   -- (i.e. \) before characters that are treated as special by
   -- pattern matcher. (The final dash is needed to avoid the
   -- regression system getting confused with other tests that have
   -- the same prefix)
   local testPrefix = string.gsub(
      string.sub(test, 3, -5), "(%W)", "%%%1") .. "%-"

   local passed, count = true, 0
   for fn in lfs.dir(outDirName) do
      local fullNm = outDirName .. fn
      local attr = lfs.attributes(fullNm)
      if attr.mode == "file" then
	 if string.find(fullNm, testPrefix) and string.sub(fullNm, -5, -1) == ".gkyl" then
	    count = count + 1 -- increment number of files compared
	    local acceptedFileNm = fullResultsDir .. "/" .. string.sub(fullNm, vloc+1, -1)
	    local status = compareFiles(acceptedFileNm, fullNm)
	    passed = passed and status
	 end
      end
   end
   -- THIS IS A HACK TO ENSURE TESTS DONT 'PASS' WHEN GKEYLL CRASHES
   if count == 0 then passed = false end
   if passed then
      numPassedTests = numPassedTests+1
      log("... passed.\n")
   else
      numFailedTests = numFailedTests+1
      log(string.format("... %s FAILED!\n", test))
   end
   return passed and 1 or 0
end

-- function to handle "run" command
local function run_action(args, name)
   loadConfigure(args)
   local luaRegTests, shellRegTests = list_tests(args)

   -- function to run after simulation is finisihed
   local postRun = function(f) return -1 end
   if args.create then
      postRun = create_action
   elseif args.check then
      postRun = check_action
   end

   log("Running regression tests ...\n\n")
   local tmStart = Time.clock()

   local function trimname(nm)
      local s,e = string.find(nm, "./")
      if s and s==1 then return string.sub(nm, e+1) end
      return nm
   end

   -- run all tests
   for _, test in pairs(luaRegTests) do
      local runtm, runlog = runLuaTest(test)
      local status = postRun(test)
      -- insert data into SQL table
      insertRegressionData(runID, trimname(test), status, runtm, runlog)
   end

   log(string.format("All regression tests completed in %g secs\n", Time.clock()-tmStart))
end

-- function to handle "listunit" command
local function listunit_action(args, name)
   loadConfigure(args)
   
   local luaUnitTests, cxxUnitTests = list_unit_tests(args)
   lume.each(luaUnitTests, print)
   lume.each(cxxUnitTests, print)
end

-- function to handle "rununit" command
local function rununit_action(args, name)
   loadConfigure(args)
   local luaUnitTests, cxxUnitTests = list_unit_tests(args)

   log("Running unit tests ...\n\n")
   local tmStart = Time.clock()
   lume.each(
      luaUnitTests, function (f) runLuaUnitTest(f) end)

   lume.each(
      cxxUnitTests, function (f) runCxxUnitTest(f) end)
   
   log(string.format("All unit tests completed in %g secs\n", Time.clock()-tmStart))
end

-- Create CLI parser to handle commands and options
local parser = argparse()
   :name("runregression")
   :require_command(true)
   :description [[
Run Gkyl regression tests. Regression tests are run and results are
compared to "accepted" results. For meaningful testing you first need
to configure the regression testing system and generate accepted
results. This can be done using the "config" and "run create"
commands. Obviously, this should be done only once.

Once accepted results are created or obtained from elsewhere,
regression testing can be  done using "run check" command.

When adding new tests it is that developer's responsibility to
ensure that the results produced are correct. Obviously, others
may not be able to determine just looking at output that the
results make sense.

The runregression tool can only be run from the gkyl/Regression
directory in the source code repo.

 ]]

parser:flag("-v --verbose", "Print verbose messages as tests are run")
parser:flag("-a --all", "Run all tests, even those in ignoretests file")

-- "config" command
local c_conf = parser:command("config", "Configure regression tests")
   :action(config_action)

c_conf:option("-p --prefix", "Location to store accepted results")
   :target("config_prefix")
c_conf:option("-m --mpiexec", "Full path to MPI executable")
   :target("config_mpiexec")
c_conf:flag("--drop-tables", "Recreate SQL tables", false)

-- "list" command
local c_list = parser:command("list", "List all regression tests")
   :action(list_action)
c_list:option("-r --run-only", "Only list this test or all tests in this directory")
c_list:option("-a --app", "Only list specified App tests. This is the prefix in test name after rt_ ")
c_list:flag("-m --moat", "Only list MOAT regression")

-- "run" tests
local c_run = parser:command("run")
   :description("Run regression tests.")
   :require_command(false)
   :action(run_action)
c_run:option("-r --run-only", "Only run these tests or all tests in these directories.\nCommma separated list")
c_run:option("-a --app", "Only run specified App tests. This is the prefix in test name after rt_ ")
c_run:flag("-m --moat", "Only run key MOAT regression")

-- check against accepted results
c_run:command("check", "Check results")

-- create accepted results
c_run:command("create", "Create accepted results")

-- "listunit" command
parser:command("listunit", "List all unit tests")
   :action(listunit_action)

-- "rununit" tests
local c_rununit = parser:command("rununit")
   :description("Run unit tests.")
   :require_command(false)
   :action(rununit_action)

-- parse command-line args (functionality encoded in command actions)
local _ = parser:parse(GKYL_COMMANDS_L)

-- print final test stats for regression tests
if numPassedTests > 0 then
   log(string.format("\nTotal tests passed: %d\n", numPassedTests))
end
if numFailedTests > 0 then
   log(string.format("Total tests FAILED: %d\n", numFailedTests))
end

if numTotalTests>0 then
   insertRegressionMeta(runID, runDate, numTotalTests, numPassedTests, numFailedTests)
end

-- ... for unit tests
if numPassedUnitTests > 0 then
   log(string.format("\nTotal unit tests passed: %d\n", numPassedUnitTests))
end
if numFailedUnitTests > 0 then
   log(string.format("Total unit tests FAILED: %d\n", numFailedUnitTests))
end
