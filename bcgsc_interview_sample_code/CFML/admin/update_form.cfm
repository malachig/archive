<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">

<html>
<head>
	<title>Record Update Form</title>
</head>

<body>

<h3>Change or Add Values and Press the Update Record Button</h3>

<!--- Retrieve the record specified by the id value passed into the template by form --->
<cfquery name="GetRecord" datasource="library_db">
	  	 SELECT entry_id, resource_type, resource_name, down_date, down_time, expected_resolution_date, expected_resolution_time, 
		 problem_class, problem, workaround, details, posted_by, date_resolved, time_resolved, post
		 FROM Lib_Db_Status
		 WHERE entry_id = #Form.entry_id#
</cfquery>

<!--- GENERATING A POP DOWN MENU WITH ALL POSSIBLE LIBRARY SUBSCRIPTION DATABASES --->
<cfquery name="GetDbName" datasource="Bibliographic databases">
		 SELECT DISTINCT title
		 FROM databases
		 ORDER BY title
</cfquery>

<!--- GENERATING A POP DOWN MENU WITH ALL POSSIBLE OTHER LIBRARY RESOURCES --->
<cfquery name="GetLibraryResource" datasource="library_db">
		 SELECT DISTINCT library_resource
		 FROM library_resources
		 ORDER BY library_resource
</cfquery>

<!--- GENERATING A POP DOWN MENU WITH ALL POSSIBLE CAMPUS RESOURCES --->
<cfquery name="GetCampusResource" datasource="library_db">
		 SELECT DISTINCT campus_resource
		 FROM campus_resources
		 ORDER BY campus_resource
</cfquery>

<!--- GENERATING A POP DOWN MENU WITH ALL POSSIBLE PROBLEM TYPES --->
<cfquery name="GetClassName" datasource="library_db">
	SELECT DISTINCT problem_class
	FROM problem_class
	ORDER BY problem_class
</cfquery>

<form action="update_record.cfm" method="post">

<!--- DATA VALIDATION --->
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

<cfoutput>
<input type="hidden" name="entry_id" value="#GetRecord.entry_id#">
<input type="hidden" name="resource_type" value="#GetRecord.resource_type#">
</cfoutput>

<table cellpadding="3" cellspacing="1">

<tr>
	<td><b>* Marked fields are required</b></td>
</tr>
<!---  
<cfoutput>
<tr>
	<td>Entry ID*</td>
	<td><input type="text" name="entry_id" size="7" maxlength="10" value="#GetRecord.entry_id#"></td>
</tr>
</cfoutput>
--->

<!--- First Choose the appropriate type of resources to display: --->
<cfif GetRecord.resource_type EQ 'Library Subscription Database'>
	<tr>
		<td>Library Subscription Database Name:*</td>
		<td><select name="resource_name">
			<cfoutput query="GetDbName">
			<option value="#title#" <cfif GetRecord.resource_name EQ GetDbName.title>SELECTED</cfif>>
			#GetDbName.title#</option>
			</cfoutput>
			</select></td>
	</tr>
</cfif>

<cfif GetRecord.resource_type EQ 'Library Resource'>
	<tr>
		<td>Library Resource Name:*</td>
		<td><select name="resource_name">
			<cfoutput query="GetLibraryResource">
			<option value="#library_resource#" <cfif GetRecord.resource_name EQ GetLibraryResource.library_resource>SELECTED</cfif>>
			#GetLibraryResource.library_resource#</option>
			</cfoutput>
			</select></td>
	</tr>
</cfif>

<cfif GetRecord.resource_type EQ 'Campus Resource'>
	<tr>
		<td>Campus Resource Name:*</td>
		<td><select name="resource_name">
			<cfoutput query="GetCampusResource">
			<option value="#campus_resource#" <cfif GetRecord.resource_name EQ GetCampusResource.campus_resource>SELECTED</cfif>>
			#GetCampusResource.campus_resource#</option>
			</cfoutput>
			</select></td>
	</tr>
</cfif>

<cfoutput>
<tr>
	<td>Down Date: (eg. 01 May 2003)*</td>
	<td><input type="text" name="down_date" size="20" maxlength="20" value="#DateFormat(GetRecord.down_date, "dd mmm yyyy")#"></td>
</tr>

<tr>
	<td>Down Time: (eg. 01:30 PM)</td>
	<td><input type="text" name="down_time" size="20" maxlength="20" value="#TimeFormat(GetRecord.down_time, "hh:mm tt")#"></td>
</tr>

<tr>
	<td>Expected Resolution Date: (eg. 01 May 2003)</td>
	<td><input type="text" name="expected_resolution_date" size="20" maxlength="20" value="#DateFormat(GetRecord.expected_resolution_date, "dd mmm yyyy")#"></td>
</tr>

<tr>
	<td>Expected Resolution Time: (eg. 01:30 PM)</td>
	<td><input type="text" name="expected_resolution_time" size="20" maxlength="20" value="#TimeFormat(GetRecord.expected_resolution_time, "hh:mm tt")#"></td>
</tr>
</cfoutput>
<tr>
	<td>Type of Downtime:*</td>
	<td><select name="problem_class">
		<cfoutput query="GetClassName">
		<option value="#problem_class#" <cfif GetRecord.problem_class EQ GetClassName.problem_class>SELECTED</cfif>>
		#GetClassName.problem_class#</option>
		</cfoutput>
		</select></td>
</tr>
<cfoutput>
<tr>
	<td>Problem:*</td>
	<td><input type="text" name="problem" size="105" maxlength="250" value="#GetRecord.problem#"></td>
</tr>

<tr>
	<td>Work Around:</td>
	<td><input type="text" name="workaround" size="105" maxlength="250" value="#GetRecord.workaround#"></td>
</tr>

<tr>
	<td>Details:</td>
	<td><input type="text" name="details" size="105" maxlength="250" value="#GetRecord.details#"></td>
</tr>

<tr>
	<td>Posted By:*</td>
	<td><input type="text" name="posted_by" size="20" maxlength="20" value="#GetRecord.posted_by#"></td>
</tr>

<tr>
	<td>Date Resolved: (eg. 01 May 2003)</td>
	<td><input type="text" name="date_resolved" size="20" maxlength="20" value="#DateFormat(GetRecord.date_resolved, "dd mmm yyyy")#"></td>
</tr>

<tr>
	<td>Time Resolved: (eg. 01:30 PM)</td>
	<td><input type="text" name="time_resolved" size="20" maxlength="20" value="#TimeFormat(GetRecord.time_resolved, "hh:mm tt")#"></td>
</tr>

<tr>
	<td>Post Record to Public Displays:*</td>
	<td><select name="post">
		<option value="Yes" <cfif GetRecord.post EQ 'Yes'>SELECTED</cfif>>Yes</option>
		<option value="No" <cfif GetRecord.post EQ 'No'>SELECTED</cfif>>No</option>
		</select></td>
</tr>
</table>
</cfoutput>

<p>Add: <input type="submit" name="update" value="Update Record"
					onclick="return confirm('Are you sure you want to update the specified record?')">
					<!--- Simple JavaScript that asks for confirmation before deleting the record --->
</form>
			

<p>
<a href="admin.cfm">Return to Admin Page</a>
<p>
<a href="..\index.cfm">Return to Database Status Page</a>
<p>
<a href="http://cybrary.uwinnipeg.ca/resources/db/index.cfm">Return to Databases</a>
</body>
</html>
