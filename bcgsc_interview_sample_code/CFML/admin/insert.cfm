<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<!--- This template inserts a record into the Lib_Db_Status table of the LibDatabases DB --->

<!--- 1.) INSERTING A LIBRARY SUBSCRIPTION DATABASE RECORD --->
<cfif IsDefined ('Form.resource_type') AND Form.resource_type EQ 'Library Subscription Database'>
<cflock name="InsertNewRecord" type="EXCLUSIVE" timeout="30">
<cftransaction>

<cfquery name="AddRecord"  datasource="library_db">
	INSERT INTO Lib_Db_Status(resource_type, resource_name, down_date, 
							  <cfif IsDate(Form.down_time)>down_time,</cfif> 
							  <cfif IsDate(Form.expected_resolution_date)>expected_resolution_date,</cfif> 
							  <cfif IsDate(Form.expected_resolution_time)>expected_resolution_time,</cfif> 
							  problem_class, problem, workaround, details, posted_by,
							  <cfif IsDate(Form.date_resolved)>date_resolved,</cfif>
							  <cfif IsDate(Form.time_resolved)>time_resolved,</cfif>
							  post)
	
	VALUES('#Form.resource_type#', '#Form.resource_name#', '#DateFormat(Form.down_date, "dd mmm yyyy")#', 
		   '#TimeFormat(Form.down_time, "hh:mm tt")#',
	       <cfif IsDate(Form.expected_resolution_date)>'#DateFormat(Form.expected_resolution_date, "dd mmm yyyy")#',</cfif> 
		   <cfif IsDate(Form.expected_resolution_time)>'#TimeFormat(Form.expected_resolution_time, "hh:mm tt")#',</cfif>
	       '#Form.problem_class#', '#Form.problem#', '#Form.workaround#', '#Form.details#', '#Form.posted_by#', 
		   <cfif IsDate(Form.date_resolved)>'#DateFormat(Form.date_resolved, "dd mmm yyyy")#',</cfif>
		   <cfif IsDate(Form.time_resolved)>'#TimeFormat(Form.time_resolved, "hh:mm tt")#',</cfif>
		   '#Form.post#')
</cfquery>
</cftransaction>
</cflock>
</cfif>


<!--- 2.) INSERTING A LIBRARY SERVER OR NETWORK RESOURCE RECORD --->
<cfif IsDefined ('Form.resource_type') AND Form.resource_type EQ 'Library Resource'>
<cflock name="InsertNewRecord" type="EXCLUSIVE" timeout="30">
<cftransaction>
<cfquery name="AddRecord"  datasource="library_db">
	INSERT INTO Lib_Db_Status(resource_type, resource_name, down_date, 
							  <cfif IsDate(Form.down_time)>down_time,</cfif> 
							  <cfif IsDate(Form.expected_resolution_date)>expected_resolution_date,</cfif> 
							  <cfif IsDate(Form.expected_resolution_time)>expected_resolution_time,</cfif> 
							  problem_class, problem, workaround, details, posted_by,
							  <cfif IsDate(Form.date_resolved)>date_resolved,</cfif>
							  <cfif IsDate(Form.time_resolved)>time_resolved,</cfif>
							  post)
	
	VALUES('#Form.resource_type#', '#Form.resource_name#', '#DateFormat(Form.down_date, "dd mmm yyyy")#', 
		   <cfif IsDate(Form.down_time)>'#TimeFormat(Form.down_time, "hh:mm tt")#',</cfif>
	       <cfif IsDate(Form.expected_resolution_date)>'#DateFormat(Form.expected_resolution_date, "dd mmm yyyy")#',</cfif> 
		   <cfif IsDate(Form.expected_resolution_time)>'#TimeFormat(Form.expected_resolution_time, "hh:mm tt")#',</cfif>
	       '#Form.problem_class#', '#Form.problem#', '#Form.workaround#', '#Form.details#', '#Form.posted_by#', 
		   <cfif IsDate(Form.date_resolved)>'#DateFormat(Form.date_resolved, "dd mmm yyyy")#',</cfif>
		   <cfif IsDate(Form.time_resolved)>'#TimeFormat(Form.time_resolved, "hh:mm tt")#',</cfif>
		   '#Form.post#')
</cfquery>
</cftransaction>
</cflock>
</cfif>


<!--- 3.) INSERTING A CAMPUS SERVER OR NETWORK RESOURCE RECORD --->
<cfif IsDefined ('Form.resource_type') AND Form.resource_type EQ 'Campus Resource'>
<cflock name="InsertNewRecord" type="EXCLUSIVE" timeout="30">
<cftransaction>
<cfquery name="AddRecord"  datasource="library_db">
	INSERT INTO Lib_Db_Status(resource_type, resource_name, down_date, 
							  <cfif IsDate(Form.down_time)>down_time,</cfif> 
							  <cfif IsDate(Form.expected_resolution_date)>expected_resolution_date,</cfif> 
							  <cfif IsDate(Form.expected_resolution_time)>expected_resolution_time,</cfif> 
							  problem_class, problem, workaround, details, posted_by,
							  <cfif IsDate(Form.date_resolved)>date_resolved,</cfif>
							  <cfif IsDate(Form.time_resolved)>time_resolved,</cfif>
							  post)
	
	VALUES('#Form.resource_type#', '#Form.resource_name#', '#DateFormat(Form.down_date, "dd mmm yyyy")#', 
		   <cfif IsDate(Form.down_time)>'#TimeFormat(Form.down_time, "hh:mm tt")#',</cfif>
	       <cfif IsDate(Form.expected_resolution_date)>'#DateFormat(Form.expected_resolution_date, "dd mmm yyyy")#',</cfif> 
		   <cfif IsDate(Form.expected_resolution_time)>'#TimeFormat(Form.expected_resolution_time, "hh:mm tt")#',</cfif>
	       '#Form.problem_class#', '#Form.problem#', '#Form.workaround#', '#Form.details#', '#Form.posted_by#', 
		   <cfif IsDate(Form.date_resolved)>'#DateFormat(Form.date_resolved, "dd mmm yyyy")#',</cfif>
		   <cfif IsDate(Form.time_resolved)>'#TimeFormat(Form.time_resolved, "hh:mm tt")#',</cfif>
		   '#Form.post#')
</cfquery>
</cftransaction>
</cflock>
</cfif>

<!--- This query retrieves the primary key value of the record we just inserted --->
<cfquery name="GetPK" datasource="library_db">
	SELECT Max(entry_id) AS MaxID
	FROM Lib_Db_Status 
</cfquery>

<!--- This query uses the value returned by the GetPK query to lookup the full record just inserted --->
<cfquery name="GetRecord" datasource="library_db">
	SELECT entry_id, resource_name, down_date, expected_resolution_date, problem, workaround, details, posted_by
	FROM Lib_Db_Status
	WHERE entry_id = #GetPK.MaxID#
</cfquery>



<html>
<head>
	<title>CFQUERY Record Insert</title>
</head>

<body>

<h2>Record Inserted!</h2>

<h3>Here are the record details ...</h3>

<table cellpadding="3" cellspacing="1">
<tr bgcolor="#888888">
	<th>Entry ID</th>
	<th>Resource Name</th>
	<th>Down Date</th>
	<th>Time</th>
	<th>Expected Resolution Date</th>
	<th>Time</th>
	<th>Problem</th>
	<th>Workaround</th>
	<th>Details</th>
	<th>Posted By</th>
</tr>

<!--- Output the new record --->
<cfoutput query="GetRecord">
<tr bgcolor="##C0C0C0">
	<td align="center">#entry_id#</td>
	<td>#resource_name#</td>
	<td align="center">#DateFormat(down_date, "dd mmm yyyy")#</td>
	<td align="center">#TimeFormat(down_time, "hh:mm tt")#</td>
	<td align="center">#DateFormat(expected_resolution_date, "dd mmm yyyy")#</td>
	<td align="center">#TimeFormat(expected_resolution_time, "hh:mm tt")#</td>
	<td>#problem#</td>
	<td>#workaround#</td>
	<td>#details#</td>
	<td align="center">#posted_by#</td>
</tr>
</cfoutput>
</table>

<p>
<h3>Add Another Record:</h3>
<form action="add_record.cfm" method="post">
<table cellpadding="3" cellspacing="1">
<tr>
	<td><input type="radio" name="Add_record" value="Library Subscription Database" CHECKED></td>
	<td>Library Subscription Database</td>
</tr>
<tr>
	<td><input type="radio" name="Add_record" value="Library Resource"></td>
	<td>Other Library Server or Network Resource</td>
</tr>
<tr>
	<td><input type="radio" name="Add_record" value="Campus Resource"></td>
	<td>Other Campus Server or Network Resource</td>
</tr>
</table>

<table cellpadding="3" cellspacing="1">
<tr>
	<td width=20></td>
	<td><input type="submit" name="Add" value="Add a Record"></td>
</tr>
</table>
</form>

<p>
<a href="admin.cfm">Return to Admin Page</a>
<p>
<a href="..\index.cfm">Return to Database Status Page</a>
<p>
<a href="http://cybrary.uwinnipeg.ca/resources/db/index.cfm">Return to Databases</a>

</body>
</html>
