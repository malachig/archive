<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">

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

<html>
<head>
	<title>Email Form</title>
</head>

<body>
<h3>Email Form - Report Database Problem Below</h3>

<pre>
Note:  If you are attempting to access library resources from off campus,
you must be a University Student or Employee.  Students must set up a proxy
service before the resources will work. 
</pre>

<form action="email_process.cfm" method="post">

<table cellpadding="3" cellspacing="1">
<tr>
	<td>Name:</td>
	<td><input type="text" name="Name" size="25" maxlength="255"></td>
</tr>

<tr>
	<td>E-Mail:</td>
	<td><input type="text" name="Email" size="25" maxlength="255"></td>
	<td>*Must be a Valid Email Address</td>
</tr>

</table>
<p>
<table cellpadding="3" cellspacing="1">
<tr>
	<td>Database Name:</td>
	<td><select name="resource_name">
		<cfoutput query="GetDbName">
		<option value="#title#" <cfif GetRecord.title EQ GetDbName.title>SELECTED</cfif>>
		<CF_SafeText2>#GetDbName.title#</CF_SafeText2></option>
		</cfoutput>
		</select></td>
</tr>
</table>

<p>
<table cellpadding="3" cellspacing="1">
<tr>
		<td><input type="Checkbox" name="OtherProblem"></td>
		<td>The resource is not listed, I will describe the problem below</td>
</tr>
</table>

<p>
Describe Problem:<br>
<textarea cols="50" rows="8" name="comments" wrap="virtual"></textarea><br>
<input type="submit" name="Submit" value="send">

</form>

<p>
<a href="index.cfm">Return to Database Status Page</a>
<p>
<a href="http://cybrary.uwinnipeg.ca/resources/db/index.cfm">Return to Databases</a>

</body>
</html>
