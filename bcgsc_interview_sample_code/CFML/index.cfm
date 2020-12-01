<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">

<!--- Get the current date and time to use as a reference --->
<cfset TheDate = Now()>
<cfset TheTime = Now()>

<html>
<head>
	<title>University of Winnipeg Library Network Resource Status Center</title>
</head>

<body>
<h2>Library Network Resource Status Center</h2>

<!--- If the user has come here with a particular database in mind.  That is if they have followed a 
	  that carried a URL variable representing the name of the resource of interest, then display 
	  information for that resource only --->
<cfif IsDefined('resourceVar')>
	<cfquery name="SpecificDb" datasource="library_db">
	SELECT resource_type, resource_name, problem, down_date, date_resolved, expected_resolution_date, expected_resolution_time, post
	FROM Lib_Db_Status
	WHERE (resource_name = '#resourceVar#' AND down_date < #TheDate#) AND ((expected_resolution_date > #TheDate#) OR (expected_resolution_date IS NULL)) 
		 		AND (resource_type = 'Library Subscription Database') AND (post = TRUE)
	ORDER BY down_date 
	</cfquery>
	
	<cfif SpecificDB.RecordCount GREATER THAN 0>	
	
	<h4>The Library Subscription Database "<cfoutput>#resourceVar#</cfoutput>" is Currently Down</h4>

	<!--- Create an HTML Table for outputting the query results.  This section
		  creates the first row of the table - used to hold the column headers.   --->
		<table cellpadding="3" cellspacing="1">
		<tr bgcolor="#888888">
			<th>DataBase Name</th>
			<th>Problem Reported</th>
			<th>Down Date</th>
			<th>Time</th>
			<th>Expected Resolution Date</th>
			<th>Time</th>
		</tr>
	
	<!--- The CFOUTPUT tag is used in conjunction with the QUERY attribute to loop
		  over each row of data in the resulting set.  During each iteration of the 
		  loop, a table row is dynamically created and populated with the query 
		  data from the current row. --->
	
		<cfoutput query="SpecificDb">
		<tr bgcolor="##C0C0C0">
			<td>#resource_name#</td>
			<td>#problem#</td>
			<td align="center">#DateFormat(down_date, "dd mmm yyyy")#</td>
			<td align="center">#TimeFormat(down_date, "hh:mm tt")#</td>
			<td align="center">#DateFormat(expected_resolution_date, "dd mmm yyyy")#</td>
			<td align="center">#TimeFormat(expected_resolution_time, "hh:mm tt")#</td>
		</tr>
		</cfoutput>
		</table>
		<cfoutput>
		<a href="output.cfm?resourceVar=#resourceVar#">View Detailed Records</a>
		</cfoutput>
		
	<cfelse>
			<h4>No Problem Currently Reported for "<cfoutput>#resourceVar#</cfoutput>"</h4>	
	</cfif>
	<p>
	<p>
	<a href="email_form.cfm">Report Problem to Library Staff</a>
	<p>
	<a href="index.cfm">Go to the Main Library Database Status Page</a>

<!--- If the user has come here via a generic link to the Status Center then display information 
	  on all databases currently known to be down and give various view options etc. --->
<cfelse>

	<!--- Query the database for Library Databases that are currently down. A database is down if the down_date
	  reported is BEFORE the current DATE (ie. the report had been entered in advance and now the time has 
	  arrived) and the date_resolved has not been entered yet (ie. it is still not functioning).  On this page
	  only records relating to Library Subscription Databases are shown.  Also they must be defined as one that
	  should be posted (ie. post=yest). 
	  *** Could set it to automatically consider the expected_resolution_date if we want to trust those times --->
	<cfquery name="GetDbName" datasource="library_db">
		 SELECT resource_type, resource_name, problem, down_date, date_resolved, expected_resolution_date, expected_resolution_time, post
		 FROM Lib_Db_Status
		 WHERE (down_date < #TheDate#) AND ((expected_resolution_date > #TheDate#) OR (expected_resolution_date IS NULL)) 
		 		AND (resource_type = 'Library Subscription Database') AND (post = TRUE)
		 ORDER BY down_date
	</cfquery>

	<!--- First query the database to get a name to display as default --->
	<cfquery name="GetRecord" datasource="Bibliographic databases">
	  	 SELECT title
		 FROM databases
		 WHERE title='Abstracts in Anthropology'
	</cfquery>

	<!--- Then query the database to get a list of all known library databases --->
	<cfquery name="GetName" datasource="Bibliographic databases">
		 SELECT DISTINCT title
		 FROM databases
		 ORDER BY title
	</cfquery>

	<h3>1. Library Subscription Databases Currently Down</h3>

	<cfif GetDbName.RecordCount GREATER THAN 0>
	<!--- Create an HTML Table for outputting the query results.  This section
		  creates the first row of the table - used to hold the column headers.   --->
		<table cellpadding="3" cellspacing="1">
		<tr bgcolor="#888888">
			<th>DataBase Name</th>
			<th>Problem Reported</th>
			<th>Down Date</th>
			<th>Time</th>
			<th>Expected Resolution Date</th>
			<th>Time</th>
		</tr>
	
	<!--- The CFOUTPUT tag is used in conjunction with the QUERY attribute to loop
		  over each row of data in the resulting set.  During each iteration of the 
		  loop, a table row is dynamically created and populated with the query 
		  data from the current row. --->
	
		<cfoutput query="GetDbName">
		<tr bgcolor="##C0C0C0">
			<td>#resource_name#</td>
			<td>#problem#</td>
			<td align="center">#DateFormat(down_date, "dd mmm yyyy")#</td>
			<td align="center">#TimeFormat(down_date, "hh:mm tt")#</td>
			<td align="center">#DateFormat(expected_resolution_date, "dd mmm yyyy")#</td>
			<td align="center">#TimeFormat(expected_resolution_time, "hh:mm tt")#</td>
		</tr>
		</cfoutput>
		</table>
	
		<table>
		<tr>
			<form action="output.cfm" method="post">
				<input type="submit" name="details" value="View Detailed Version">
			</form>
		</tr>
		</table>
			
	<cfelse>
			<h4>No Databases Currently Reported as Down</h4>
	</cfif>

	<h3>2. View History for a Specific Resource</h3>
	<form action="output.cfm" method="post">
		<table cellpadding="3" cellspacing="1">
		<tr>
			<td>Database Name:</td>
		</tr>
		<tr>
			<td><select name="resource_name">
				<cfoutput query="GetName">
				<option value="#title#" <cfif GetRecord.title EQ GetName.title>SELECTED</cfif>>
				<CF_SafeText2>#GetName.title#</CF_SafeText2></option>
				</cfoutput>
				</select></td>
		</tr>
		</table>
		<input type="submit" name="single" value="View Records">
	</form>

	<form action="output.cfm" method="post">

	<!--- Radio Button Section, allowing user to specify if they want more detailed information about
		  websites that are currently down.  Choices:
	   	  1.) View history for last 2 weeks 
	 	  2.) View history for last month
	 	  3.) View reports affecting next 2 weeks
	  	  4.) View reports affecting next month
	 	  5.) View detailed report of those currently down, instead of the summary seen above --->

	<h3>3. Advanced Queries</h3>

	<b>Make a Selection and Press Submit Button:</b>
	<table width=700>

		<tr>
			<td colspan="2"><b>Choose the Desired Time Frame</b></td>
			<td colspan="2"><b>Choose the Type of Resource</b></td>
		</tr>
		<tr>
			<td align="right"><input type="Radio" name="Selection" value="next_2_weeks" CHECKED></td>
			<td>Next 2 Weeks</td>
			<td><input type="Radio" name="Resource" value="Library Subscription Database" CHECKED></td>
			<td>Library Subscription Database</td>
		</tr>
		<tr>
			<td align="right"><input type="Radio" name="Selection" value="next_month"></td>
			<td>Next Month</td>
			<td><input type="Radio" name="Resource" value="Library Resource"></td>
			<td>Other Library Server or Network Resources</td>
		</tr>
		<tr>
			<td align="right"><input type="Radio" name="Selection" value="last_2_weeks"></td>
			<td>Last 2 Weeks</td>
			<td><input type="Radio" name="Resource" value="Campus Resource"></td>
			<td>Other Campus Server or Network Resources</td>
		</tr>
		<tr>
			<td align="right"><input type="Radio" name="Selection" value="last_month"></td>
			<td>Last Month</td>
	</table>
	<input type="submit" name="submit" value="Submit Query">
	</form>

	<h3>4. Instructions:</h3>
<pre>
If you are having trouble accessing a particular database or library resource 
check the list above for the databases that are currently known to be down.

If you would like more information about databases or library resources that may be down 
or will be down in the future, make a selection above and press the 'submit query' button.
	
If a database is not working and you can find no reference to the problem, please click the 
"report problem" link to report this problem to library staff.  Include the name of the 
database and the nature of problem.
</pre>

<h3>5. Administration:</h3>

	<!--- Link to Admin page so that new reports can be made and database status can be updated --->
	<a href="admin\authenticate.cfm">Go to Admin Pages</a>
	
	<p>
	<a href="email_form.cfm">Report Problem to Library Staff</a>

	
	
</cfif>

<p>
<a href="http://cybrary.uwinnipeg.ca/resources/db/index.cfm">Return to Databases</a>

</body>
</html>
