-- Gkyl ------------------------------------------------------------------------
--
-- Tool to handle gkyl help system.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local argparse = require "Lib.argparse"
local keywordMap = require "Tool.HelpKeywordMap"

-- Create CLI parser to handle commands and options
local parser = argparse()
   :name("man")
   :description [[
Open Gkeyll documentation webpage at gkeyll.rtfd.io and optionally
search within it via keywords. To just open the documentation page
type:

 gkyl man

For example, to open documentation for basis functions do:

 gkyl man basis

or

 gkyl man modal

The keywords don't need to be completely spelled out. For example the
following will also work:

 gkyl man bas

If there are more than one keywork matches the tool will display
options for the user to choose from. If no keyword matches are found
the home page is opened.
]]

parser:flag("-k --keywords", "List available keywords")
parser:option("-s --search", "Search documentation instead")

parser:argument("optional", "Keyword to search")
   :args("?")

-- parse command line parameters
local args = parser:parse(GKYL_COMMANDS_L)

if args.keywords then
   kwds = {}
   for i, kv in ipairs(keywordMap) do
      kwds[i] = kv[1]
   end
   table.sort(kwds)
   for _, kv in ipairs(kwds) do
      io.write(string.format("'%s' ", kv))
   end
   
   print("")
   return
end

-- Attempts to open a given URL in the system default browser
local open_cmd
local function open_url(url)
   if not open_cmd then
      if package.config:sub(1,1) == '\\' then -- windows
	 open_cmd = function(url)
	    -- Should work on anything since (and including) win'95
	    os.execute(string.format('start "%s"', url))
	 end
      elseif io.popen("uname -s"):read'*a' == "Darwin\n" then -- OSX/Darwin
	 open_cmd = function(url)
	    os.execute(string.format('open "%s"', url))
	 end
      else
	 open_cmd = function(url)
	    os.execute(string.format('xdg-open "%s"', url))
	 end
      end
   end
   open_cmd(url)
end

-- list of matching keywords
local function getMatchingKeywords(keyword)
   local matches = {}
   local count = 1
   for _, kv in ipairs(keywordMap) do
      if string.match(kv[1], keyword) then
	 matches[count] = kv
	 count = count+1
      end
   end
   return matches
end

if not args.optional then
   open_url("https://gkeyll.readthedocs.io/en/latest/")
else
   local kw = args.optional
   local matches = getMatchingKeywords(kw:lower())

   -- take action based on how many matches were returned
   if #matches == 0 then
      print(string.format("No such keyword '%s'. Launching search.", args.optional))
      open_url(string.format(
		  "https://gkeyll.readthedocs.io/en/latest/search.html?q=%s", args.optional)
      )
   elseif #matches == 1 then
      open_url(matches[1][2])
   else
      -- display matches to user
      print("Multiple matches. Choose one:")
      for i, ma in ipairs(matches) do
	 print(string.format("%d: %s", i, ma[1]))
      end
      io.write(">: ")
      local n = io.read("*number")
      if n <= #matches then
	 -- open corresponding URL
	 open_url(matches[n][2])
      end
   end
end
