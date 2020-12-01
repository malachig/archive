<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">

<!--- *** Note: Change this to email address or adresses that you wont problems reported to. *** --->
<cfset ContactEmail = "mgriffit@io.uwinnipeg.ca">

<!--- CASE 1: If a problem was reported for a resources OTHER than one of the library subscription DBs 
	  listed.  ie. the checkbox to that effect was checked. --->
<cfif IsDefined ('Form.OtherProblem')>
<cfmail FROM="#Form.Email#"
	    TO=#ContactEmail#
		SUBJECT="Library Resource Problem Report">
The following information was submitted via a feedback form:

Name: #Form.Name#
E-mail: #Form.Email#
Comments:
#Form.comments#
</cfmail>

<html>
<head>
	<title>Feedback Form Processor</title>
</head>

<body>

<h2>The resource problem has been reported successfully!</h2>
The following information was submitted:<br>
<cfoutput>
Name: #Form.Name#<br>
E-Mail: #Form.Email#<br>
Comments:<br>
#Form.comments#
</cfoutput>



<!--- If a problem was reported for one of the library subscription DBs listed.  ie. the checkbox
	  to that effect was NOT checked. --->
<cfelse>
<cfmail FROM="#Form.Email#"
	    TO=#ContactEmail#
		SUBJECT="Database Problem Report for #Form.resource_name#">
The following information was submitted via a feedback form:

Name: #Form.Name#
E-mail: #Form.Email#
Database: #Form.resource_name#
Comments:
#Form.comments#
</cfmail>

<html>
<head>
	<title>Feedback Form Processor</title>
</head>

<body>

<h2>The database problem has been reported successfully!</h2>
The following information was submitted:<br>
<cfoutput>
Name: #Form.Name#<br>
E-Mail: #Form.Email#<br>
Database: #Form.resource_name#<br>
Comments:<br>
#Form.comments#
</cfoutput>
</cfif>

<p>
<a href="index.cfm">Return to Database Status Page</a>
<p>
<a href="http://cybrary.uwinnipeg.ca/resources/db/index.cfm">Return to Databases</a>

</body>
</html>
