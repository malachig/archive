<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">

<html>
<head>
	<title>Creating URL Parameters</title>
</head>

<body>

<!--- Set variables to be passed to another template --->
<cfset  resourceVar="Abstracts in Anthropology">

<h2>Passing Data via URL Parameters</h2>

<!--- Create a hyperlink containing URL parameters --->

<cfoutput>
<a href="index.cfm?resourceVar=#resourceVar#">Click this link to pass the URL parameters</a>
</cfoutput>

</body>
</html>
