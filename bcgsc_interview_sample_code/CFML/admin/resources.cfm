<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">

<!--- Get current list of library resources --->
<cfquery name="LibraryResources" datasource="library_db">
		 SELECT ID, library_resource
		 FROM library_resources
		 ORDER BY library_resource 
</cfquery>

<!--- Get current list of campus resources --->
<cfquery name="CampusResources" datasource="library_db">
	   	 SELECT ID, campus_resource
		 FROM campus_resources
		 ORDER BY campus_resource
</cfquery>

<html>
<head>
	<title>Managing Library and Campus Resources</title>
</head>

<body>



<!--- ADDING RESOURCES
	  Selecting the type of resource, filling in its name and submitting the record
	  to the database. --->
<cfif IsDefined ('Form.resources') AND Form.resources EQ 'Add Resource'>

<h2>Adding a Library or Campus Resource</h2>

<pre>
This option allows you to add new resources for tracking purposes.  After a new resource is
added it will appear in the list of resources available for reporting an outage.
</pre>

<!--- The following section deals with adding LIBRARY or CAMPUS resources to the list of those that
may be tracked.  First show those currently available (Grouped by type).  Then allow the user to fill 
in a form (select type and fill in its name) and press the add button.  
Then go to the insert_record.cfm page. --->

<h3>Current Library Resources</h3>
<cfif LibraryResources.RecordCount GREATER THAN 0>
	<table cellpadding="3" cellspacing="1">
	<tr bgcolor="#888888">
		<th>Library Resource Name</th>
	</tr>
	<cfoutput query="LibraryResources">
	<tr bgcolor="##C0C0C0">
		<td>#library_resource#</td>
	</tr>
	</cfoutput>
	</table>
<cfelse>
	<h4>No Library Resources have Been Entered</h4>
</cfif>

<h3>Current Campus Resources</h3>
<cfif CampusResources.RecordCount GREATER THAN 0>
	<table cellpadding="3" cellspacing="1">
	<tr bgcolor="#888888">
		<th>Campus Resource Name</th>
	</tr>
	<cfoutput query="CampusResources">
	<tr bgcolor="##C0C0C0">
	<td>#campus_resource#</td>
	</tr>
	</cfoutput>
	</table>
<cfelse>
	<h4>No Campus Resources have Been Entered</h4>
</cfif>
<p>
Choose the Type of Resource, enter the resource name and press the 'Add Resource' Button

<form action="insert_resource.cfm" method="post">

<!--- DATA VALIDATION --->
<input type="hidden" name="resource_name_Required" value="Resource name is a required field">
<input type="hidden" name="resource_type_Required" value="Resource type is a required field">

<table cellpadding="3" cellspacing="1">

<tr>
	<td><b>* Marked fields are required</b></td>
</tr>

<tr>
	<td>Resource Type:*</td>
	<td><select name="resource_type">
		<option value="Library Resource" SELECTED>Library Resource</option>
		<option value="Campus Resource">Campus Resource</option>
		</select>
	</td>
</tr>

<tr>
	<td>Resource Name:*</td>
	<td><input type="text" name="resource_name" size="40" maxlength="100"></td>
</tr>
</table>

<p><input type="submit" name="submit" value="Add Resource"
onclick="return confirm('Are you sure you want to add this record?')">
	
</form>
</cfif>

<!--- DELETING RESOURCES
	  Selecting the resource by name and submitting the delete request. --->
<cfif (IsDefined('Form.resources') AND (Form.resources EQ 'Delete Resource')) 
	   OR (IsDefined('Form.resources') AND (Form.resources EQ 'Delete Library Resource'))
	   OR (IsDefined('Form.resources') AND (Form.resources EQ 'Delete Campus Resource'))>

<cfif IsDefined('Form.resources') AND Form.resources EQ 'Delete Library Resource'>
	<cfquery NAME="DeleteRecord" datasource="library_db">
         DELETE FROM library_resources
		 WHERE ID = #ID#
	</cfquery>
</cfif>

<cfif IsDefined('Form.resources') AND Form.resources EQ 'Delete Campus Resource'>
	<cfquery NAME="DeleteRecord" datasource="library_db">
         DELETE FROM campus_resources
		 WHERE ID = #ID#
	</cfquery>
</cfif>

<!--- Redo Query to Get An Update for Populating the selection boxes. --->
<!--- Get current list of library resources --->
<cfquery name="LibraryResources" datasource="library_db">
		 SELECT ID, library_resource
		 FROM library_resources
		 ORDER BY library_resource 
</cfquery>

<!--- Get current list of campus resources --->
<cfquery name="CampusResources" datasource="library_db">
	   	 SELECT ID, campus_resource
		 FROM campus_resources
		 ORDER BY campus_resource
</cfquery>

<h2>Deleting a Library or Campus Resource</h2>

<h3>Library Resources</h3>

<cfif LibraryResources.RecordCount GREATER THAN 0>
	<form action="resources.cfm" method="post">
	<table border="0">
	<tr>
		<td>
		<select name="ID" size="5">
			<cfoutput query="LibraryResources">
				<option value="#ID#">#library_resource#</option>
			</cfoutput>
		</select>
		</td>
	</tr>

	<tr>
		<td align="center">
			<input type="submit" name="resources" value="Delete Library Resource" 
				onclick="return confirm('Are you sure you want to delete the specified record?')">
				<!--- Simple JavaScript that asks for confirmation before deleting the record --->
		</td>
	</tr>
	</table>
	</form>
<cfelse>
	<h4>No Library Resources have Been Entered</h4>
</cfif>

<h3>Campus Resources</h3>

<cfif CampusResources.RecordCount GREATER THAN 0>
	<form action="resources.cfm" method="post">
	<table border="0">
	<tr>
		<td>
		<select name="ID" size="5">
			<cfoutput query="CampusResources">
				<option value="#ID#">#campus_resource#</option>
			</cfoutput>
		</select>
		</td>
	</tr>
	
	<tr>
		<td align="center">
			<input type="submit" name="resources" value="Delete Campus Resource" 
				onclick="return confirm('Are you sure you want to delete the specified record?')">
				<!--- Simple JavaScript that asks for confirmation before deleting the record --->
		</td>
	</tr>
	</table>
	</form>
<cfelse>
	<h4>No Campus Resources have Been Entered</h4>
</cfif>

</cfif>
<p>
<a href="admin.cfm">Return to Admin Page</a>
<p>
<a href="..\index.cfm">Return to Database Status Page</a>
<p>
<a href="http://cybrary.uwinnipeg.ca/resources/db/index.cfm">Return to Databases</a>

</body>
</html>
