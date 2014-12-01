<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<!--- This page will be used by administrators to update the database as reports from library
	  databases are received.  PASSWORD AUTHENTICATION will be needed to limit who can access THIS PAGE --->
<html>
<head>
	<title>University of Winnipeg Library Database Administration Page</title>
</head>

<body>

<h2>Adminstration Page</h2>

<H3>Select Which Action You Wish to Perform:</H3>

<h3>Add Records</h3>
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

<h3>Update Records</h3>
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
</form>

<h3>View Records</h3>
<form action="view_all.cfm" method="post">
<table cellpadding="3" cellspacing="1">
<tr>
	<td width=20></td>
	<td><input type="submit" name="ViewAll" value="View All Records"></td>
</tr>
</table>
</form>

<h3>Add or Delete Library and Campus Resources </h3>
<form action="resources.cfm" method="post">
<table cellpadding="3" cellspacing="1">
<tr>
	<td><input type="radio" name="resources" value="Add Resource" CHECKED></td>
	<td>Add Other Library/Campus Server or Network Resources</td>
</tr>
<tr>
	<td><input type="radio" name="resources" value="Delete Resource"></td>
	<td>Delete Other Library/Campus Server or Network Resources</td>
</tr>
</table>

<table cellpadding="3" cellspacing="1">
<tr>
	<td width=20></td>
	<td width><input type="submit" name="submit" value="Add or Delete Resources"></td>
</tr>
</table>
</form>

<h3>Delete Records</h3>
<form action="delete_record.cfm" method="post">
<table cellpadding="3" cellspacing="1">
<tr>
	<td width=20></td>
	<td><input type="submit" name="Remove" value="Delete a Record"></td>
</tr>
</table>
</form>

<p>
<a href="..\index.cfm">Return to Database Status Page</a>
<p>
<a href="http://cybrary.uwinnipeg.ca/resources/db/index.cfm">Return to Databases</a>
</body>
</html>
