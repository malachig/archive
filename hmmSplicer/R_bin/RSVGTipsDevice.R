# load library
library("RSVGTipsDevice")

# Define the SVG characteristics
devSVGTips(svg_filename, toolTipMode=2, title=svg_title, toolTipFontSize=10,sub.special=FALSE,height=8,width=12)

# Create sensitive areas
# Mousover text

setSVGShapeToolTip(title=title_text,desc1=first_line_text,desc2=second_line_text)
# you can ignore desc1 and desc2 and just have text # desc2 is only available when  toolTipMode=2 in SVG declaration


# Url
# target="_blank" makes the URL open in a new tab setSVGShapeURL(url_text, target="_blank")

# plot the element
# Ex1 text:
        mtext(Title, outer = FALSE, cex = 1)

# Ex2 an invisible area
         rect(i-0.4,y.max,i+0.4,y.min,col=NULL,density=NULL,border=NA)




# Warning, if you declare the clickable area and then you plot something over it, that new plot will not be clickable.
# For example, if you make a rectangle clickable and then plot a normal triangle inside, the area of the rectangle overlapping the triangle will not be clickable