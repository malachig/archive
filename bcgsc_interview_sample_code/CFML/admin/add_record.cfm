<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">

<!--- THIS PAGE DOES ONE OF THREE THINGS DEPENDING ON WHICH TYPE OF RECORD IS BEING ADDED, AS SELECTED
	  ON THE CALLING PAGE "admin.cfm" 
	  Note: although this makes this page very large and repetitive, doing it this way will allow me to 
	  make the three types of entries different later.  For example we may want to require different fields
	  for the different types of records, etc.--->

<!--- Get the current date and time to use as a reference --->
<cfset TheDate = Now()>
<cfset TheTime = Now()>
	  
<html>
<head>
	<title>Adding a Record to the Database</title>
</head>

<body>

<!--- GENERATING POP DOWN MENU WITH DIFFERENT PROBLEM CLASSES --->
<cfquery name="GetClass" datasource="library_db">
	SELECT problem_class
	FROM problem_class
	WHERE problem_class = 'Scheduled Downtime'
</cfquery>

<cfquery name="GetClassName" datasource="library_db">
	SELECT DISTINCT problem_class
	FROM problem_class
	ORDER BY problem_class
</cfquery>


<!--- 1.)  ADDING A LIBRARY SUBSCRIPTION DATABASE RECORD--->

<cfif IsDefined ('Form.Add_record') AND Form.Add_record EQ 'Library Subscription Database'>
<!--- GENERATING A POP DOWN MENU WITH ALL POSSIBLE LIBRARY SUBSCRIPTION DATABASES --->

<!--- First query the database to get a name to display as default --->
<cfquery name="GetRecord" datasource="Bibliographic databases">
	  	 SELECT title
		 FROM databases
		 WHERE title='Abstracts in Anthropology'
</cfquery>

<!--- Then query the database to get a list of all known library databases --->
<cfquery name="GetDbName" datasource="Bibliographic databases">
		 SELECT DISTINCT title
		 FROM databases
		 ORDER BY title
</cfquery>

<h3>Adding a Library Subscription Database Record</h3>

<form action="insert.cfm" method="post">

<!--- DATA VALIDATION --->
<!--- Note: There is a problem with data checking for dates.  If you leave the date or time blank
	  this causes an error because a blank field is not a valid date or time.  To avoid this problem, in the 
	  insert.cfm file I use a series of IsDate() expressions to include the date in the insert query if it is 
	  defined and not if it was left blank.  This way users can leave the field blank but if they do specify a 
	  date or time the format is checked for errors. --->

<input type="hidden" name="resource_name_Required" value="Resource name is a required field">
<input type="hidden" name="down_date_Date" value="Down Date must be a valid date in date format">
<input type="hidden" name="down_date_Required" value="Down Date is Required">
<input type="hidden" name="down_time_Time" value="Down time must be in time format">
<input type="hidden" name="expected_resolution_date_Date" value="Expected resolution date must be a valid date in date format">
<input type="hidden" name="expected_resolution_time_Time" value="Expected resolution time must be in time format">
<input type="hidden" name="problem_class_Required" value="Type of problem must be specified">
<input type="hidden" name="problem_Required" value="Nature of problem is a required field">
<input type="hidden" name="date_resolved_Date" value="Date Resolved must be a valid date in date format">
<input type="hidden" name="time_resolved_Time" value="Time Resolved must be in time format">
<input type="hidden" name="posted_by_Required" value="Posted by is a required field">

<table cellpadding="3" cellspacing="1">

<tr>
	<td><b>* Marked fields are required</b></td>
</tr>

<tr>
	<td>Database Name:*</td>
	<td><select name="resource_name">
		<cfoutput query="GetDbName">
		<option value="#title#" <cfif GetRecord.title EQ GetDbName.title>SELECTED</cfif>>
		<CF_SafeText2>#GetDbName.title#</CF_SafeText2></option>
		</cfoutput>
		</select></td>
</tr>
<tr>
	<td>Down Date: (eg. 01 May 2003)*</td>
	<td><input type="text" name="down_date" size="20" maxlength="20" value=<cfoutput>"#DateFormat(TheDate, "dd mmm yyyy")#"</cfoutput>></td>
</tr>

<tr>
	<td>Down Time: (eg. 01:30 PM)</td>
	<td><input type="text" name="down_time" size="20" maxlength="20" value=<cfoutput>"#TimeFormat(TheTime, "hh:mm tt")#"</cfoutput>></td>
</tr>

<tr>
	<td>Expected Resolution Date: (eg. 01 May 2003)</td>
	<td><input type="text" name="expected_resolution_date" size="20" maxlength="20"></td>
</tr>

<tr>
	<td>Expected Resolution Time: (eg. 01:30 PM)</td>
	<td><input type="text" name="expected_resolution_time" size="20" maxlength="20"></td>
</tr>

<tr>
	<td>Type of Downtime:*</td>
	<td><select name="problem_class">
		<cfoutput query="GetClassName">
		<option value="#problem_class#" <cfif GetClass.problem_class EQ GetClassName.problem_class>SELECTED</cfif>>
		#GetClassName.problem_class#</option>
		</cfoutput>
		</select></td>
</tr>

<tr>
	<td>Problem:*</td>
	<td><input type="text" name="problem" size="105" maxlength="250"></td>
</tr>

<tr>
	<td>Work Around:</td>
	<td><input type="text" name="workaround" size="105" maxlength="250"></td>
</tr>

<tr>
	<td>Details:</td>
	<td><input type="text" name="details" size="105" maxlength="250"></td>
</tr>

<tr>
	<td>Posted By:*</td>
	<td><input type="text" name="posted_by" size="20" maxlength="20"></td>
</tr>

<tr>
	<td>Date Resolved: (eg. 01 May 2003)</td>
	<td><input type="text" name="date_resolved" size="20" maxlength="20"></td>
</tr>

<tr>
	<td>Time Resolved: (eg. 01:30 PM)</td>
	<td><input type="text" name="time_resolved" size="20" maxlength="20"></td>
</tr>

<tr>
	<td>Post Record to Public Displays:*</td>
	<td><select name="post">
		<option value="Yes" SELECTED>Yes</option>
		<option value="No">No</option>
		</select></td>
</tr>
</table>

<p>Add Record: <input type="submit" name="resource_type" value="Library Subscription Database"
				onclick="return confirm('Are you sure the information above is correct?')">
</form>
</cfif>


<!--- 2.)  ADDING AN OTHER LIBRARY RESOURCE RECORD --->

<cfif IsDefined ('Form.Add_record') AND Form.Add_record EQ 'Library Resource'>
<!--- GENERATING A POP DOWN MENU WITH ALL POSSIBLE OTHER LIBRARY RESOURCES --->

<cfquery name="GetRecord" datasource="library_db">
	  	 SELECT library_resource
		 FROM library_resources
		 WHERE library_resource='Cybrary'
</cfquery>
<cfquery name="GetResourceName" datasource="library_db">
		 SELECT DISTINCT library_resource
		 FROM library_resources
		 ORDER BY library_resource
</cfquery>

<h3>Adding a Library Server or Network Resource Record</h3>

<form action="insert.cfm" method="post">

<input type="hidden" name="resource_name_Required" value="Resource name is a required field">
<input type="hidden" name="down_date_Date" value="Down Date must be a valid date in date format">
<input type="hidden" name="down_date_Required" value="Down Date is Required">
<input type="hidden" name="down_time_Time" value="Down time must be in time format">
<input type="hidden" name="expected_resolution_date_Date" value="Expected resolution date must be a valid date in date format">
<input type="hidden" name="expected_resolution_time_Time" value="Expected resolution time must be in time format">
<input type="hidden" name="problem_class_Required" value="Type of problem must be specified">
<input type="hidden" name="problem_Required" value="Nature of problem is a required field">
<input type="hidden" name="date_resolved_Date" value="Date Resolved must be a valid date in date format">
<input type="hidden" name="time_resolved_Time" value="Time Resolved must be in time format">
<input type="hidden" name="posted_by_Required" value="Posted by is a required field">

<table cellpadding="3" cellspacing="1">

<tr>
	<td><b>* Marked fields are required</b></td>
</tr>

<tr>
	<td>Library Resource Name:*</td>
	<td><select name="resource_name">
		<cfoutput query="GetResourceName">
		<option value="#library_resource#" <cfif GetRecord.library_resource EQ GetResourceName.library_resource>SELECTED</cfif>>
		#GetResourceName.library_resource#</option>
		</cfoutput>
		</select></td>
</tr>
<tr>
	<td>Down Date: (eg. 01 May 2003)*</td>
	<td><input type="text" name="down_date" size="20" maxlength="20" value=<cfoutput>"#DateFormat(TheDate, "dd mmm yyyy")#"</cfoutput>></td>
</tr>

<tr>
	<td>Down Time: (eg. 01:30 PM)</td>
	<td><input type="text" name="down_time" size="20" maxlength="20" value=<cfoutput>"#TimeFormat(TheTime, "hh:mm tt")#"</cfoutput>></td>
</tr>

<tr>
	<td>Expected Resolution Date: (eg. 01 May 2003)</td>
	<td><input type="text" name="expected_resolution_date" size="20" maxlength="20"></td>
</tr>

<tr>
	<td>Expected Resolution Time: (eg. 01:30 PM)</td>
	<td><input type="text" name="expected_resolution_time" size="20" maxlength="20"></td>
</tr>

<tr>
	<td>Type of Downtime:*</td>
	<td><select name="problem_class">
		<cfoutput query="GetClassName">
		<option value="#problem_class#" <cfif GetClass.problem_class EQ GetClassName.problem_class>SELECTED</cfif>>
		#GetClassName.problem_class#</option>
		</cfoutput>
		</select></td>
</tr>

<tr>
	<td>Problem:*</td>
	<td><input type="text" name="problem" size="105" maxlength="250"></td>
</tr>

<tr>
	<td>Work Around:</td>
	<td><input type="text" name="workaround" size="105" maxlength="250"></td>
</tr>

<tr>
	<td>Details:</td>
	<td><input type="text" name="details" size="105" maxlength="250"></td>
</tr>

<tr>
	<td>Posted By:*</td>
	<td><input type="text" name="posted_by" size="20" maxlength="20"></td>
</tr>

<tr>
	<td>Date Resolved: (eg. 01 May 2003)</td>
	<td><input type="text" name="date_resolved" size="20" maxlength="20"></td>
</tr>

<tr>
	<td>Time Resolved: (eg. 01:30 PM)</td>
	<td><input type="text" name="time_resolved" size="20" maxlength="20"></td>
</tr>

<tr>
	<td>Post Record to Public Displays:*</td>
	<td><select name="post">
		<option value="Yes" SELECTED>Yes</option>
		<option value="No">No</option>
		</select></td>
</tr>
</table>

<p>Add Record: <input type="submit" name="resource_type" value="Library Resource"
				onclick="return confirm('Are you sure the information above is correct?')">
</form>
</cfif>

<!--- 3.)  ADDING A CAMPUS RESOURCE RECORD --->

<cfif IsDefined ('Form.Add_record') AND Form.Add_record EQ 'Campus Resource'>

<!--- GENERATING A POP DOWN MENU WITH ALL POSSIBLE CAMPUS RESOURCES --->
<cfquery name="GetRecord" datasource="library_db">
	  	 SELECT campus_resource
		 FROM campus_resources
		 WHERE campus_resource='IO Server'
</cfquery>
<cfquery name="GetResourceName" datasource="library_db">
		 SELECT DISTINCT campus_resource
		 FROM campus_resources
		 ORDER BY campus_resource
</cfquery>

<h3>Adding a Campus Server or Network Resource Record</h3>

<form action="insert.cfm" method="post">

<input type="hidden" name="resource_name_Required" value="Resource name is a required field">
<input type="hidden" name="down_date_Date" value="Down Date must be a valid date in date format">
<input type="hidden" name="down_date_Required" value="Down Date is Required">
<input type="hidden" name="down_time_Time" value="Down time must be in time format">
<input type="hidden" name="expected_resolution_date_Date" value="Expected resolution date must be a valid date in date format">
<input type="hidden" name="expected_resolution_time_Time" value="Expected resolution time must be in time format">
<input type="hidden" name="problem_class_Required" value="Type of problem must be specified">
<input type="hidden" name="problem_Required" value="Nature of problem is a required field">
<input type="hidden" name="date_resolved_Date" value="Date Resolved must be a valid date in date format">
<input type="hidden" name="time_resolved_Time" value="Time Resolved must be in time format">
<input type="hidden" name="posted_by_Required" value="Posted by is a required field">

<table cellpadding="3" cellspacing="1">

<tr>
	<td><b>* Marked fields are required</b></td>
</tr>

<tr>
	<td>Campus Resource Name:*</td>
	<td><select name="resource_name">
		<cfoutput query="GetResourceName">
		<option value="#campus_resource#" <cfif GetRecord.campus_resource EQ GetResourceName.campus_resource>SELECTED</cfif>>
		#GetResourceName.campus_resource#</option>
		</cfoutput>
		</select></td>
</tr>
<tr>
	<td>Down Date: (eg. 01 May 2003)*</td>
	<td><input type="text" name="down_date" size="20" maxlength="20" value=<cfoutput>"#DateFormat(TheDate, "dd mmm yyyy")#"</cfoutput>></td>
</tr>

<tr>
	<td>Down Time: (eg. 01:30 PM)</td>
	<td><input type="text" name="down_time" size="20" maxlength="20" value=<cfoutput>"#TimeFormat(TheTime, "hh:mm tt")#"</cfoutput>></td>
</tr>

<tr>
	<td>Expected Resolution Date: (eg. 01 May 2003)</td>
	<td><input type="text" name="expected_resolution_date" size="20" maxlength="20"></td>
</tr>

<tr>
	<td>Expected Resolution Time: (eg. 01:30 PM)</td>
	<td><input type="text" name="expected_resolution_time" size="20" maxlength="20"></td>
</tr>

<tr>
	<td>Type of Downtime:*</td>
	<td><select name="problem_class">
		<cfoutput query="GetClassName">
		<option value="#problem_class#" <cfif GetClass.problem_class EQ GetClassName.problem_class>SELECTED</cfif>>
		#GetClassName.problem_class#</option>
		</cfoutput>
		</select></td>
</tr>

<tr>
	<td>Problem:*</td>
	<td><input type="text" name="problem" size="105" maxlength="250"></td>
</tr>

<tr>
	<td>Work Around:</td>
	<td><input type="text" name="workaround" size="105" maxlength="250"></td>
</tr>

<tr>
	<td>Details:</td>
	<td><input type="text" name="details" size="105" maxlength="250"></td>
</tr>

<tr>
	<td>Posted By:*</td>
	<td><input type="text" name="posted_by" size="20" maxlength="20"></td>
</tr>

<tr>
	<td>Date Resolved: (eg. 01 May 2003)</td>
	<td><input type="text" name="date_resolved" size="20" maxlength="20"></td>
</tr>

<tr>
	<td>Time Resolved: (eg. 01:30 PM)</td>
	<td><input type="text" name="time_resolved" size="20" maxlength="20"></td>
</tr>

<tr>
	<td>Post Record to Public Displays:*</td>
	<td><select name="post">
		<option value="Yes" SELECTED>Yes</option>
		<option value="No">No</option>
		</select></td>
</tr>
</table>

<p>Add Record: <input type="submit" name="resource_type" value="Campus Resource"
				onclick="return confirm('Are you sure the information above is correct?')">
</form>
</cfif>
<p>
<a href="admin.cfm">Return to Admin Page</a>
<p>
<a href="..\index.cfm">Return to Database Status Page</a>
<p>
<a href="http://cybrary.uwinnipeg.ca/resources/db/index.cfm">Return to Databases</a>

</body>
</html>
