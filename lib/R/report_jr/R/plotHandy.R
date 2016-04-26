
plotHandy <- function(p, legtxt, legcols, legtitle) {

# legtxt <- as.character(levels(assay1$gentreatment))
# legcols <- c("violet", "orange", "gray40", "black")
# legtitle <- "group"

fracPlot <- 0.9 # fraction of total plot area allocated to the main plot

# Define 2 viewports
vp1x <- fracPlot/2
fudge <- 0.5 # needed for aesthetics to get legend centered in white space
vp2x <- fudge*((1-fracPlot)/2) + fracPlot
vp1 <- viewport(x = vp1x, y = 0.5, width = fracPlot, height = 1.0) # plot space
vp2 <- viewport(x = vp2x, y = 0.5, width = 1 - fracPlot, height = 1.0) # legend space
	
# Create a legend

leg <- draw.key(key = list(text = list(lab = legtxt,  columns = 1, col = legcols), cex.title = 1.0,
	title = legtitle))

# Now plot it all

grid.newpage()
pushViewport(vp1)
print(p, newpage = FALSE)
popViewport()
pushViewport(vp2)
grid.draw(leg)

}
