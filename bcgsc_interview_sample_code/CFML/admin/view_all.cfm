<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">

<!--- This query retrieves the complete record of all entries in the database --->
<cfif IsDefined('Form.SortSelection')>
	<cfquery name="ViewAll" datasource="library_db">
		SELECT entry_id, resource_name, down_date, down_time, expected_resolution_date, expected_resolution_time, problem, 
	    	   workaround, details, posted_by, date_resolved, time_resolved
		FROM Lib_Db_Status
		ORDER BY #Form.sort#
	</cfquery>

<cfelse>
	<cfquery name="ViewAll" datasource="library_db">
		SELECT entry_id, resource_name, down_date, down_time, expected_resolution_date, expected_resolution_time, problem, 
	    	   workaround, details, posted_by, date_resolved, time_resolved
		FROM Lib_Db_Status
		ORDER BY entry_id
	</cfquery>
</cfif>


<html>
<head>
	<title>Complete Listing of All Records in the ScoreBoard Database</title>
</head>

<body>

<cfif ViewAll.RecordCount GREATER THAN 0>
	<h3>Full record details ...</h3>

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
		<th>Date Resolved</th>
		<th>Time</th>
	</tr>

	<!--- Output the new record --->
	<cfoutput query="ViewAll">
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
		<td align="center">#DateFormat(date_resolved, "dd mmm yyyy")#</td>
		<td align="center">#TimeFormat(time_resolved, "hh:mm tt")#</td>
	</tr>
	</cfoutput>
	</table>
	
	<form action="view_all.cfm" method="post">
	<table cellpadding="3" cellspacing="1">
	<tr>
		Sort Records By:
	</tr>
	<tr>
		<td><select name="sort">
			<option value="entry_id" SELECTED>Entry ID</option>
			<option value="resource_name">Name</option>
			<option value="down_date">Down Date</option>
			<option value="expected_resolution_date">Expected Resolution Date</option>
			<option value="posted_by">Posted By</option>
			</select></td>
	</tr>
	</table>
	
	<input type="submit" name="SortSelection" value="Sort">
	</form>
<cfelse>
	<h4>The Database is Empty</h4>
</cfif>

<a href="admin.cfm">Return to Admin Page</a>
<p>
<a href="..\index.cfm">Return to Database Status Page</a>
<p>
<a href="http://cybrary.uwinnipeg.ca/resources/db/index.cfm">Return to Databases</a>

</body>
</html>
