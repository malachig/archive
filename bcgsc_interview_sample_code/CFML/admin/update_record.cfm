<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">

<!--- Update the record specifed in the entry_id field.  Note that numeric 
	  values are not enclosed in single quotes in the set clause. --->

<cfquery name="UpdateRecord" Datasource="library_db">
		 UPDATE Lib_Db_Status
		 SET resource_type = '#Form.resource_type#',
    		 resource_name = '#Form.resource_name#',
			 <cfif IsDate(Form.down_date)>down_date = '#DateFormat(Form.down_date, "dd mmm yyyy")#',</cfif>
			 <cfif IsDate(Form.down_time)>down_time = '#TimeFormat(Form.down_time, "hh:mm tt")#',</cfif>
			 <cfif IsDate(Form.expected_resolution_date)>expected_resolution_date = '#DateFormat(Form.expected_resolution_date, "dd mmm yyyy")#',</cfif>
			 <cfif IsDate(Form.expected_resolution_time)>expected_resolution_time = '#TimeFormat(Form.expected_resolution_time, "hh:mm tt")#',</cfif>
			 problem_class = '#Form.problem_class#',
			 problem = '#Form.problem#',
			 workaround = '#Form.workaround#',
			 details = '#Form.details#',
			 posted_by = '#Form.posted_by#',
			 <cfif IsDate(Form.date_resolved)>date_resolved = '#DateFormat(Form.date_resolved, "dd mmm yyyy")#',</cfif>
			 <cfif IsDate(Form.time_resolved)>time_resolved = '#TimeFormat(Form.time_resolved, "hh:mm tt")#',</cfif>
			 post = '#Form.post#'
		 WHERE entry_id = #Form.entry_id#
</cfquery>


<!--- Retrieve the record that was just updated --->
<cfquery name="GetRecord" datasource="library_db">
	     SELECT entry_id, resource_type, resource_name, down_date, down_time, expected_resolution_date,
		        expected_resolution_time, problem_class, problem, workaround, details, posted_by,
				date_resolved, time_resolved, post
		 FROM Lib_Db_Status
		 WHERE entry_id = #Form.entry_id#
</cfquery>

<html>
<head>
	<title>Record Update</title>
</head>

<body>

<h2>Record Updated!</h2>

<h3>Here are the record details ...</h3>

<table cellpadding="3" cellspacing="1">
<tr bgcolor="#888888">
	<th>Entry ID</th>
	<th>Resource Type</th>
	<th>Resource Name</th>
	<th>Down Date</th>
	<th>Time</th>
	<th>Expected Resolution Date</th>
	<th>Time</th>
	<th>Problem Type</th>
	<th>Problem</th>
	<th>Workaround</th>
	<th>Details</th>
	<th>Posted By</th>
	<th>Date Resolved</th>
	<th>Time</th>
	<th>Post</th>
</tr>

<!--- Output the record --->
<cfoutput query="GetRecord">
<tr bgcolor="##C0C0C0">
	<td>#entry_id#</td>
	<td>#resource_type#</td>
	<td>#resource_name#</td>
	<td>#DateFormat(down_date, "dd mmm yyyy")#</td>
	<td>#TimeFormat(down_time, "hh:mm tt")#</td>
	<td>#DateFormat(expected_resolution_date, "dd mmm yyyy")#</td>
	<td>#TimeFormat(expected_resolution_time, "hh:mm tt")#</td>
	<td>#problem_class#</td>
	<td>#problem#</td>
	<td>#workaround#</td>
	<td>#details#</td>
	<td>#posted_by#</td>
	<td>#DateFormat(date_resolved, "dd mmm yyyy")#</td>
	<td>#TimeFormat(time_resolved, "hh:mm tt")#</td>
	<td>#post#</td>
	
</tr>
</cfoutput>
</table>

<h3>Update Another Records</h3>
<form action="choose_update.cfm" method="post">
<table cellpadding="3" cellspacing="1">
<tr>
	<td><input type="radio" name="Update_record" value="Show All" CHECKED></td>
	<td>Show All Entries</td>
</tr>
<tr>
	<td><input type="radio" name="Update_record" value="Library Subscription Database"></td>
	<td>Library Subscription Databases Only</td>
</tr>
<tr>
	<td><input type="radio" name="Update_record" value="Library Resource"></td>
	<td>Other Library Server or Network Resources</td>
</tr>
<tr>
	<td><input type="radio" name="Update_record" value="Campus Resource"></td>
	<td>Other Campus Server or Network Resources</td>
</tr>
</table>

<table cellpadding="3" cellspacing="1">
<tr>
	<td width=20></td>
	<td width><input type="submit" name="Update" value="Update a Record"></td>
</tr>
</table>


<p>
<a href="admin.cfm">Return to Admin Page</a>
<p>
<a href="..\index.cfm">Return to Database Status Page</a>
<p>
<a href="http://cybrary.uwinnipeg.ca/resources/db/index.cfm">Return to Databases</a>
</body>
</html>
