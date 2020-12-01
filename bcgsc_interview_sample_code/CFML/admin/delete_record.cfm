<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">

<!--- Check to see if the form is calling itself to delete a single record.
	  Note that by NOT scoping the ID variable in the WHERE clause, we Allow
	  the template to accept both Form and URL variables. --->
<cfif IsDefined('Delete')>
	<cfquery NAME="DeleteRecord" datasource="library_db">
         DELETE FROM Lib_Db_Status
		 WHERE entry_id = #entry_id#
	</cfquery>
</cfif>

<!--- This query retrieves basic information about all records that have not been resolved
	  so they can be displayed. --->
	  
<cfset TheDate = Now()>
<cfquery name="DisplayList" datasource="library_db">
	SELECT entry_id, resource_name, problem, posted_by, down_date
	FROM Lib_Db_Status
	WHERE date_resolved IS NULL
	ORDER BY down_date
</cfquery>

<html>
<head>
	<title>Deleting a Record from the Database</title>
</head>

<body>

<h2>Please Select the Record to Delete from the Database</h2>
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
	<h3>Only records that have an unresolved status are displayed.</h3>
		Choose the record to delete from the list and press "Delete"
	<p>
	<form action="delete_record.cfm" method="post">
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
			<input type="submit" name="delete" value="delete" 
				onclick="return confirm('Are you sure you want to delete the specified record?')">
				<!--- Simple JavaScript that asks for confirmation before deleting the record --->
		</td>
	</tr>
	</table>
	</form>
<cfelse>
		<h4>No Records to be Deleted</h4>
		<pre>
		Note: Records that have a resolution date associated with them, ie. the status of 
		the report has been resolved, are not eligable for deletion.
		</pre>
</cfif>

<p>
<a href="admin.cfm">Return to Admin Page</a>
<p>
<a href="..\index.cfm">Return to Database Status Page</a>
<p>
<a href="http://cybrary.uwinnipeg.ca/resources/db/index.cfm">Return to Databases</a>

</body>
</html>
