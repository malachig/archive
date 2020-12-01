<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">

<cfset shadow="claudius">

<html>
<head>
	<title>Administration Password Authentication Page</title>
</head>

<body>

<Form action="authenticate.cfm" method="post">
<table cellpadding="3" cellspacing="1">
<tr>
	<td>Password:</td>
	<td><input type="text" name="passwd" size="20" maxlength="25"</td>
</tr>

<tr>
	<td colspan="2"><input type="submit" name="Submit" value="Submit"</td>
</tr>

</table>
</form>

<cfif IsDefined('Form.passwd') AND Form.passwd EQ shadow>
	Password Accepted --->
	<a href="admin.cfm">You May Continue</a>
<cfelseif IsDefined('Form.passwd')>
	Password Incorrect ...
</cfif>

<p>
<a href="..\index.cfm">Return to Database Status Page</a>
<p>
<a href="http://cybrary.uwinnipeg.ca/resources/db/index.cfm">Return to Databases</a>

</body>
</html>
