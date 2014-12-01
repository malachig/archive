<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">

<cfset TheDate = Now()>

<cfif IsDefined('Form.Update_record') AND Form.Update_record EQ 'Show All'>
	<cfquery name="DisplayList" datasource="library_db">
		SELECT entry_id, resource_name, problem, posted_by, down_date
		FROM Lib_Db_Status
		ORDER BY resource_name
	</cfquery>
</cfif>

<cfif IsDefined('Form.Update_record') AND Form.Update_record EQ 'Library Subscription Database'>
	<cfquery name="DisplayList" datasource="library_db">
	SELECT entry_id, resource_type, resource_name, problem, posted_by, down_date
	FROM Lib_Db_Status
	WHERE resource_type = 'Library Subscription Database'
	ORDER BY resource_name
	</cfquery>
</cfif>

<cfif IsDefined('Form.Update_record') AND Form.Update_record EQ 'Library Resource'>
	<cfquery name="DisplayList" datasource="library_db">
	SELECT entry_id, resource_type, resource_name, problem, posted_by, down_date
	FROM Lib_Db_Status
	WHERE resource_type = 'Library Resource'
	ORDER BY resource_name
	</cfquery>
</cfif>

<cfif IsDefined('Form.Update_record') AND Form.Update_record EQ 'Campus Resource'>
	<cfquery name="DisplayList" datasource="library_db">
	SELECT entry_id, resource_type, resource_name, problem, posted_by, down_date
	FROM Lib_Db_Status
	WHERE resource_type = 'Campus Resource'
	ORDER BY resource_name
	</cfquery>
</cfif>


<html>
<head>
	<title>Choosing a Record to Update</title>
</head>

<body>

<h2>Please Select the Record to Update from the Database</h2>

<cfif DisplayList.RecordCount GREATER THAN 0>
	<table cellpadding="3" cellspacing="1">
	<tr bgcolor="#888888">
		<th>Record ID</th>
		<th align="left">DataBase Name</th>
		<th>Problem Reported</th>
		<th>Posted By</th>
		<th>Down Date</th>
	</tr>

	<!--- The CFOUTPUT tag is used in conjunction with the QUERY attribute to loop
		  over each row of data in the resulting set.  During each iteration of the 
		  loop, a table row is dynamically created and populated with the query 
		  data from the current row. --->
	
	<cfoutput query="DisplayList">
	<tr bgcolor="##C0C0C0">
		<td align="center">#entry_id#</td>
		<td>#resource_name#</td>
		<td>#problem#</td>
		<td align="center">#posted_by#</td>
		<td aligh="center">#DateFormat(down_date, "dd mmm yyyy")#</td>
	</tr>
	</cfoutput>
	</table>

	<p>
	<form action="update_form.cfm" method="post">
	<table border="0">
	<tr>
		<td>
		<select name="entry_id" size="5">
			<cfoutput query="DisplayList">
				<option value="#entry_id#">#entry_id# -> #resource_name#</option>
			</cfoutput>
		</select>
		</td>
	</tr>

	<tr>
		<td align="center">
			<input type="submit" name="update" value="update">
		</td>
	</tr>
	</table>
	</form>

<cfelse>
	<h4>No Records Meet the Selected Criteria</h4>
</cfif>

<p>
<a href="admin.cfm">Return to Admin Page</a>
<p>
<a href="..\index.cfm">Return to Database Status Page</a>
<p>
<a href="http://cybrary.uwinnipeg.ca/resources/db/index.cfm">Return to Databases</a>
</body>
</html>
