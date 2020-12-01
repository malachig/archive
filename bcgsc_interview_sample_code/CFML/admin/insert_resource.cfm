<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">

<cfif IsDefined ('Form.resource_type') AND Form.resource_type EQ 'Library Resource'>
<cflock name="InsertNewResource" type="EXCLUSIVE" timeout="30">
<cftransaction>
<cfquery name="AddResource"  datasource="library_db">
	INSERT INTO library_resources(library_resource)
	VALUES('#Form.resource_name#')
</cfquery>
</cftransaction>
</cflock>
</cfif>

<cfif IsDefined ('Form.resource_type') AND Form.resource_type EQ 'Campus Resource'>
<cflock name="InsertNewResource" type="EXCLUSIVE" timeout="30">
<cftransaction>
<cfquery name="AddResource"  datasource="library_db">
	INSERT INTO campus_resources(campus_resource)
	VALUES('#Form.resource_name#')
</cfquery>
</cftransaction>
</cflock>
</cfif>

<!--- Retrieves the primary key value of the record we just inserted --->
<cfif IsDefined ('Form.resource_type') AND Form.resource_type EQ 'Library Resource'>
<cfquery name="GetPK" datasource="library_db">
	SELECT Max(ID) AS MaxID
	FROM library_resources
</cfquery>
<!--- This query uses the value returned by the GetPK query to lookup the full record just inserted --->
<cfquery name="GetRecord" datasource="library_db">
	SELECT ID, library_resource
	FROM library_resources
	WHERE ID = #GetPK.MaxID#
</cfquery>
</cfif>

<cfif IsDefined ('Form.resource_type') AND Form.resource_type EQ 'Campus Resource'>
<cfquery name="GetPK" datasource="library_db">
	SELECT Max(ID) AS MaxID
	FROM campus_resources
</cfquery>
<cfquery name="GetRecord" datasource="library_db">
	SELECT ID, campus_resource
	FROM campus_resources
	WHERE ID = #GetPK.MaxID#
</cfquery>
</cfif>

<html>
<head>
	<title>Resource Insert</title>
</head>

<body>

<h2>Record Inserted!</h2>

<h3>Here are the record details ...</h3>

<cfif IsDefined ('Form.resource_type') AND Form.resource_type EQ 'Library Resource'>
<table cellpadding="3" cellspacing="1">
<tr bgcolor="#888888">
	<th>ID</th>
	<th>Library Resource Name</th>
</tr>

<!--- Output the new record --->
<cfoutput query="GetRecord">
<tr bgcolor="##C0C0C0">
	<td align="center">#ID#</td>
	<td>#library_resource#</td>
</tr>
</cfoutput>
</table>
</cfif>

<cfif IsDefined ('Form.resource_type') AND Form.resource_type EQ 'Campus Resource'>
<table cellpadding="3" cellspacing="1">
<tr bgcolor="#888888">
	<th>ID</th>
	<th>Campus Resource Name</th>
</tr>

<!--- Output the new record --->
<cfoutput query="GetRecord">
<tr bgcolor="##C0C0C0">
	<td align="center">#ID#</td>
	<td>#campus_resource#</td>
</tr>
</cfoutput>
</table>
</cfif>

<h3>Add Another Resource:</h3>
<form action="resources.cfm" method="post">
<table cellpadding="3" cellspacing="1">
<tr>
	<td><input type="submit" name="resources" value="Add Resource"></td>
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
