<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">

<!--- Set up date variables needed for various queries --->
<cfset TheDate = Now()>
<cfset TwoWeeksAgo = DateAdd('ww', -2, TheDate)>
<cfset MonthAgo = DateAdd('m', -1, TheDate)>
<cfset TwoWeeksAhead = DateAdd('ww', 2, TheDate)>
<cfset MonthAhead = DateAdd('m', 1, TheDate)>

<!--- What happens on this page will be determined by what what was selected on the db_query page.
	  For performance reasons, each query will be contained within CFIF statements so that the
	  database is only accessed as much as neccessary. --->

<html>
<head>
	<title>Resource Status Report Page</title>
</head>

<body>

<!--- If the "Detailed Records" selection was made, display the appropriate query results --->	
<cfif IsDefined ('Form.details')>
	<cfquery name="detailed_records" datasource="library_db">
	SELECT entry_id, resource_name, down_date, down_time, expected_resolution_date, expected_resolution_time, problem, 
	       workaround, details, posted_by, date_resolved, time_resolved
	FROM Lib_Db_Status
	WHERE (down_date < #TheDate#) AND ((expected_resolution_date > #TheDate#) OR (expected_resolution_date IS NULL)) 
		   AND (resource_type = 'Library Subscription Database') AND (post = TRUE)
	ORDER BY down_date
	</cfquery>
<h3>Full record details ...</h3>

	<cfif detailed_records.RecordCount GREATER THAN 0>
		<table cellpadding="3" cellspacing="1">
		<tr bgcolor="#888888">
			<th>Entry ID</th>
			<th>Library Subscription Database</th>
			<th>Down Date</th>
			<th>Time</th>
			<th>Expected Resolution Date</th>
			<th>Time</th>
			<th>Problem</th>
			<th>Workaround</th>
			<th>Details</th>
			<th>Posted By</th>
			<th>Date Resolved</th>
			<th>Time</th>
		</tr>
		
		<!--- Output the new record --->
		<cfoutput query="detailed_records">
		<tr bgcolor="##C0C0C0">
			<td align="center">#entry_id#</td>
			<td>#resource_name#</td>
			<td align="center">#DateFormat(down_date, "dd mmm yyyy")#</td>
			<td align="center">#TimeFormat(down_time, "hh:mm tt")#</td>
			<td align="center">#DateFormat(expected_resolution_date, "dd mmm yyyy")#</td>
			<td align="center">#TimeFormat(expected_resolution_time, "hh:mm tt")#</td>
			<td>#problem#</td>
			<td>#workaround#</td>
			<td>#details#</td>
			<td align="center">#posted_by#</td>
			<td align="center">#DateFormat(date_resolved, "dd mmm yyyy")#</td>
			<td align="center">#TimeFormat(time_resolved, "hh:mm tt")#</td>
		</tr>
		</cfoutput>
		</table>
	<cfelse>
		<h4>No Records Meet the Selected Criteria</h4>
	</cfif>
</cfif>

<!--- If the "View History for a Specific Resource" selection was made, display the appropriate query results --->	
<cfif IsDefined ('Form.single')>
	<cfquery name="single_resource" datasource="library_db">
	SELECT entry_id, resource_name, down_date, down_time, expected_resolution_date, expected_resolution_time, problem, 
	       workaround, details, posted_by, date_resolved, time_resolved
	FROM Lib_Db_Status
	WHERE resource_name = '#Form.resource_name#'
	ORDER BY down_date
	</cfquery>
<h3>Full record details ...</h3>

	<cfif single_resource.RecordCount GREATER THAN 0>
		<table cellpadding="3" cellspacing="1">
		<tr bgcolor="#888888">
			<th>Entry ID</th>
			<th>Library Subscription Database</th>
			<th>Down Date</th>
			<th>Time</th>
			<th>Expected Resolution Date</th>
			<th>Time</th>
			<th>Problem</th>
			<th>Workaround</th>
			<th>Details</th>
			<th>Posted By</th>
			<th>Date Resolved</th>
			<th>Time</th>
		</tr>
		
		<!--- Output the new record --->
		<cfoutput query="single_resource">
		<tr bgcolor="##C0C0C0">
			<td align="center">#entry_id#</td>
			<td>#resource_name#</td>
			<td align="center">#DateFormat(down_date, "dd mmm yyyy")#</td>
			<td align="center">#TimeFormat(down_time, "hh:mm tt")#</td>
			<td align="center">#DateFormat(expected_resolution_date, "dd mmm yyyy")#</td>
			<td align="center">#TimeFormat(expected_resolution_time, "hh:mm tt")#</td>
			<td>#problem#</td>
			<td>#workaround#</td>
			<td>#details#</td>
			<td align="center">#posted_by#</td>
			<td align="center">#DateFormat(date_resolved, "dd mmm yyyy")#</td>
			<td align="center">#TimeFormat(time_resolved, "hh:mm tt")#</td>
		</tr>
		</cfoutput>
		</table>
	<cfelse>
		<h4>No Records Meet the Selected Criteria</h4>
	</cfif>
</cfif>

	
	
<!--- If the "Last 2 Weeks" selection was made, display the appropriate query results 
      Either it went down during this period or was still down from before the previous 
	  two weeks --->	
<cfif IsDefined ('Form.Selection') AND Form.Selection EQ 'last_2_weeks'>
	<cfquery name="last_2_weeks" datasource="library_db">
	SELECT resource_name, problem, down_date, date_resolved, expected_resolution_date, expected_resolution_time
	FROM Lib_Db_Status
	WHERE (((down_date < #TheDate#) AND (down_date > #TwoWeeksAgo#)) OR ((date_resolved > #TwoWeeksAgo#) AND (date_resolved < #TheDate#)))
	       AND (resource_type = '#Form.Resource#') AND (post = TRUE)
	ORDER BY down_date
	</cfquery>

	
	<cfoutput><h3>#Form.Resource#s Down in the Last Two Weeks</h3></cfoutput>
	<cfif last_2_weeks.RecordCount GREATER THAN 0>
		<table cellpadding="3" cellspacing="1">
		<tr bgcolor="#888888">
			<th><cfoutput>#Form.Resource#</cfoutput></th>
			<th>Problem Reported</th>
			<th>Down Date</th>
			<th>Time</th>
			<th>Expected Resolution Date</th>
			<th>Time</th>
		</tr>

		<cfoutput query="last_2_weeks">	
		<tr bgcolor="##C0C0C0">
			<td>#resource_name#</td>
			<td>#problem#</td>
			<td align="center">#DateFormat(down_date, "dd mmm yyyy")#</td>
			<td align="center">#TimeFormat(down_date, "hh:mm tt")#</td>
			<td align="center">#DateFormat(expected_resolution_date, "dd mmm yyyy")#</td>
			<td align="center">#TimeFormat(expected_resolution_time, "hh:mm tt")#</td>
		</tr>
		</cfoutput>
		</table>
	<cfelse>
		<h4>No Records Meet the Selected Criteria</h4>
	</cfif>
</cfif>


<!--- If the "Last Month" selection was made, display the appropriate query results --->	
<cfif IsDefined ('Form.Selection') AND Form.Selection EQ 'last_month'>
	<cfquery name="last_month" datasource="library_db">
	SELECT resource_name, problem, down_date, date_resolved, expected_resolution_date, expected_resolution_time
	FROM Lib_Db_Status
	WHERE (((down_date < #TheDate#) AND (down_date > #MonthAgo#)) OR ((date_resolved > #MonthAgo#) AND (date_resolved < #TheDate#)))
		   AND (resource_type = '#Form.Resource#') AND (post = TRUE)
	ORDER BY down_date
	</cfquery>

	<cfoutput><h3>#Form.Resource#s Down in the Last Month</h3></cfoutput>
	<cfif last_month.RecordCount GREATER THAN 0>
		<table cellpadding="3" cellspacing="1">
		<tr bgcolor="#888888">
			<th><cfoutput>#Form.Resource#</cfoutput></th>
			<th>Problem Reported</th>
			<th>Down Date</th>
			<th>Time</th>
			<th>Expected Resolution Date</th>
			<th>Time</th>
		</tr>
		<cfoutput query="last_month">
		<tr bgcolor="##C0C0C0">
			<td>#resource_name#</td>
			<td>#problem#</td>
			<td align="center">#DateFormat(down_date, "dd mmm yyyy")#</td>
			<td align="center">#TimeFormat(down_date, "hh:mm tt")#</td>
			<td align="center">#DateFormat(expected_resolution_date, "dd mmm yyyy")#</td>
			<td align="center">#TimeFormat(expected_resolution_time, "hh:mm tt")#</td>
		</tr>
		</cfoutput>
		</table>
	<cfelse>
		<h4>No Records Meet the Selected Criteria</h4>
	</cfif>
</cfif>


<!--- If the "Next 2 Weeks" selection was made, display the appropriate query results --->	
<cfif IsDefined ('Form.Selection') AND Form.Selection EQ 'next_2_weeks'>
	<cfquery name="next_2_weeks" datasource="library_db">
	SELECT resource_name, problem, down_date, date_resolved, expected_resolution_date, expected_resolution_time
	FROM Lib_Db_Status
	WHERE (down_date < #TwoWeeksAhead#) AND (expected_resolution_date > #TheDate#) AND (resource_type = '#Form.Resource#') AND (post = TRUE)
	ORDER BY down_date
	</cfquery>

	<cfoutput><h3>#Form.Resource#s Expected to Have Problems in the Next Two Weeks</h3></cfoutput>
	<cfif next_2_weeks.RecordCount GREATER THAN 0>
		<table cellpadding="3" cellspacing="1">
		<tr bgcolor="#888888">
			<th><cfoutput>#Form.Resource#</cfoutput></th>
			<th>Problem Reported</th>
			<th>Down Date</th>
			<th>Time</th>
			<th>Expected Resolution Date</th>
			<th>Time</th>
		</tr>
		<cfoutput query="next_2_weeks">
		<tr bgcolor="##C0C0C0">
			<td>#resource_name#</td>
			<td>#problem#</td>
			<td align="center">#DateFormat(down_date, "dd mmm yyyy")#</td>
			<td align="center">#TimeFormat(down_date, "hh:mm tt")#</td>
			<td align="center">#DateFormat(expected_resolution_date, "dd mmm yyyy")#</td>
			<td align="center">#TimeFormat(expected_resolution_time, "hh:mm tt")#</td>
		</tr>
		</cfoutput>
		</table>
	<cfelse>
		<h4>No Records Meet the Selected Criteria</h4>
	</cfif>
</cfif>

<!--- If the "Next Month" selection was made, display the appropriate query results --->	
<cfif IsDefined ('Form.Selection') AND Form.Selection EQ 'next_month'>
	<cfquery name="next_month" datasource="library_db">
	SELECT resource_name, problem, down_date, date_resolved, expected_resolution_date, expected_resolution_time
	FROM Lib_Db_Status
	WHERE (down_date < #MonthAhead#) AND (expected_resolution_date > #TheDate#) AND (resource_type = '#Form.Resource#') AND (post = TRUE)
	ORDER BY down_date
	</cfquery>

	<cfoutput><h3>#Form.Resource#s Expected to Have Problems in the Next Month</h3></cfoutput>
	<cfif next_month.RecordCount GREATER THAN 0>
		<table cellpadding="3" cellspacing="1">
		<tr bgcolor="#888888">
			<th><cfoutput>#Form.Resource#</cfoutput></th>
			<th>Problem Reported</th>
			<th>Down Date</th>
			<th>Time</th>
			<th>Expected Resolution Date</th>
			<th>Time</th>
		</tr>
		<cfoutput query="next_month">
		<tr bgcolor="##C0C0C0">
			<td>#resource_name#</td>
			<td>#problem#</td>
			<td align="center">#DateFormat(down_date, "dd mmm yyyy")#</td>
			<td align="center">#TimeFormat(down_date, "hh:mm tt")#</td>
			<td align="center">#DateFormat(expected_resolution_date, "dd mmm yyyy")#</td>
			<td align="center">#TimeFormat(expected_resolution_time, "hh:mm tt")#</td>
		</tr>
		</cfoutput>
		</table>
	<cfelse>
		<h4>No Records Meet the Selected Criteria</h4>
	</cfif>
</cfif>

<!--- If the user entered the status page carrying a variable for the name of a single resource.  
	  It is forwared to this page for the detailed display  --->

<cfif IsDefined ('resourceVar')>
	<cfquery name="SpecificDb_detailed" datasource="library_db">
	SELECT entry_id, resource_name, down_date, down_time, expected_resolution_date, expected_resolution_time, problem, 
	       workaround, details, posted_by, date_resolved, time_resolved
	FROM Lib_Db_Status
	WHERE (resource_name = '#resourceVar#') AND (resource_type = 'Library Subscription Database') AND (post = TRUE)
	ORDER BY down_date
	</cfquery>
<h3>Full record details ...</h3>

	<cfif SpecificDb_detailed.RecordCount GREATER THAN 0>
		<table cellpadding="3" cellspacing="1">
		<tr bgcolor="#888888">
			<th>Entry ID</th>
			<th>Library Subscription Database</th>
			<th>Down Date</th>
			<th>Time</th>
			<th>Expected Resolution Date</th>
			<th>Time</th>
			<th>Problem</th>
			<th>Workaround</th>
			<th>Details</th>
			<th>Posted By</th>
			<th>Date Resolved</th>
			<th>Time</th>
		</tr>
		
		<!--- Output the new record --->
		<cfoutput query="SpecificDb_detailed">
		<tr bgcolor="##C0C0C0">
			<td align="center">#entry_id#</td>
			<td>#resource_name#</td>
			<td align="center">#DateFormat(down_date, "dd mmm yyyy")#</td>
			<td align="center">#TimeFormat(down_time, "hh:mm tt")#</td>
			<td align="center">#DateFormat(expected_resolution_date, "dd mmm yyyy")#</td>
			<td align="center">#TimeFormat(expected_resolution_time, "hh:mm tt")#</td>
			<td>#problem#</td>
			<td>#workaround#</td>
			<td>#details#</td>
			<td align="center">#posted_by#</td>
			<td align="center">#DateFormat(date_resolved, "dd mmm yyyy")#</td>
			<td align="center">#TimeFormat(time_resolved, "hh:mm tt")#</td>
		</tr>
		</cfoutput>
		</table>
	<cfelse>
		<h4>No Records Meet the Selected Criteria</h4>
	</cfif>
</cfif>


<p>
<a href="index.cfm">Return to Database Status Page</a>
<p>
<a href="http://cybrary.uwinnipeg.ca/resources/db/index.cfm">Return to Databases</a>

</body>
</html>
