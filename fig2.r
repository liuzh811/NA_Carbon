#plot together
d1.ec = data.frame(read.csv("F:/zhihua/dataset/results2/EC.sensitivity.csv"))[,-1]
d1.rs = data.frame(read.csv("F:/zhihua/dataset/results2/RS.sensitivity.csv"))[,-1]
d1.trendy = data.frame(read.csv("F:/zhihua/dataset/results2/trendy.sensitivity.csv"))[,-1]
d1.lpj = data.frame(read.csv("F:/zhihua/dataset/results2/lpj.sensitivity.csv"))[,-1]

D1 = rbind(data.frame(d1.ec, Method = "EC"),
     data.frame(d1.rs, Method = "RS"),
     data.frame(d1.trendy, Method = "TRENDY"),
     data.frame(d1.trendy, Method = "LPJ"))
 
D1 = rbind(data.frame(d1.ec, Method = "EC"),
     data.frame(d1.rs, Method = "RS"),
     data.frame(d1.trendy, Method = "TRENDY"))
 

levels(D1$type)  <- c("Spatial","Temporal")
levels(D1$region)  <- c("< 750","> 750")

ylab=expression("" ~ g ~ C ~ m^{-2} ~ yr ^{-1}~ "/100mm")

p1 <- ggplot(data=D1, aes(x=region, y=mean, fill=flux)) + 
    geom_bar(colour="black", stat="identity",
             position=position_dodge(),
             size=.3) + 
    facet_grid(type~Method)+			 
	ylab(ylab) + 	 
    # ylab(expression(paste(beta ["temporal"]))) + 
	# xlab("Region by MAP") + 
	# ylab("") + 
	xlab("MAP (mm)") + 
    # ggtitle("Average bill for 2 people") +     # Set title
    theme_bw() + 
	geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) + 
	theme(legend.position=c(.85, .4)) + 	
	theme(legend.title=element_blank()) +
	theme(legend.text = element_text(size = 18)) +
	theme(axis.title.x = element_text(face="bold", colour="black", size=18),axis.text.x  = element_text(colour="black",size=18))+
    theme(axis.title.y = element_text(face="bold", colour="black", size=18),axis.text.y  = element_text(colour="black",size=18))+
    theme(strip.text.x = element_text(size=18))+
    theme(strip.text.y = element_text(size=18)) 

color1 = c("grey20", "grey90")
p1 + scale_fill_manual(values=color1, 
                    name="",
                    breaks=levels(D1$flux),
                    labels=levels(D1$flux)) 

ggsave("F:/zhihua/dataset/results2/sensitivity.bw.png", width = 8, height = 6, units = "in")

color2 = c("blue", "red")
p1 + scale_fill_manual(values=color2, 
                    name="",
                    breaks=levels(D1$flux),
                    labels=levels(D1$flux)) 

ggsave("F:/zhihua/dataset/results2/sensitivity.color.png", width = 8, height = 6, units = "in")

####new figure 2
#plot together
d.temporal.ec = data.frame(read.csv("F:/zhihua/dataset/results2/delt.temporal.ec.csv"))[,-1]
d.spatial.ec = data.frame(read.csv("F:/zhihua/dataset/results2/delt.spatial.ec.csv"))[,-1]

d.temporal.rs = data.frame(read.csv("F:/zhihua/dataset/results2/delt.temporal.rs.csv"))[,-1]
d.spatial.rs = data.frame(read.csv("F:/zhihua/dataset/results2/delt.spatial.rs.csv"))[,-1]

d.temporal.trendy = data.frame(read.csv("F:/zhihua/dataset/results2/delt.temporal.trendy.csv"))[,-1]
d.spatial.trendy = data.frame(read.csv("F:/zhihua/dataset/results2/delt.spatial.trendy.csv"))[,-1]

prep.grd = seq(100,1300, length.out = 25)

png("F:/zhihua/dataset/results2/fig2-4.png",height = 2500, width = 2000, res = 300, units = "px")

par(mfrow=c(2,1),mar=c(0,0,0,0),oma=c(5,5,0,0))

 # plot temporal sensitivity
# http://htmlcolorcodes.com/
# par(mar=c(0,0,0,0))

plot(1, type="n", 
								  ylim = c(-100,100),
                                  				  xlim = c(100, 1300),								  
								  cex = 2, lwd = 4, 
								  col = "green",
								  xlab = "Anuual Mean Precipitation (MAP: mm)", 
								  # ylab = "Spatial sensitivity", 
								  # ylab = expression(paste(beta ["temporal"])), 
								  cex.axis = 1.5, 
								  cex.lab = 1.6,
								  xaxt='n', ann=FALSE,yaxt='n')

# add axis and labels
axis(side = 2, las = 1, tck = -.02, at = c(-50,0,50), labels = c("-50","0","50"), cex.axis = 1.5)

abline(h = 0, lty = 2, lwd = 2)		
								  
d.temporal.ec = d.temporal.ec[order(d.temporal.ec$prep),]
mylm = 	lm(mean~prep, data = d.temporal.ec)
abline(mylm,col="green", lwd = 3)

prd<-predict(mylm,newdata=data.frame(prep = prep.grd),interval = c("confidence"), level = 0.50,type="response") 
polygon(c(rev(prep.grd), prep.grd), 
		c(rev(prd[,2]), prd[,3]), 
        col=rgb(0, 0.5, 0,0.25),
		border = NA)
	
abline(v = 885, col = "green", lwd = 4, lty = 2)
# abline(v = 885-85, col = "green")
# abline(v = 885+85, col = "green")

x1 = c(885-85,885+85)
y1 = c(500,500)
y2 = c(-500,-500)
# abline(v = c(650, 750, 700 + 100), lty = c(2,1,2))
# abline(v = c(650, 700 + 100), lty = c(2,2))

polygon(c(x1, rev(x1)), c(y1, rev(y2)),
      col=rgb(0, 0.5, 0,0.25), 
border = NA)


# overlay remote sensing	
d.temporal.rs = d.temporal.rs[order(d.temporal.rs$prep),]
mylm = 	lm(mean~prep, data = d.temporal.rs)
abline(mylm,col="blue", lwd = 3)

prd<-predict(mylm,newdata=data.frame(prep = prep.grd),interval = c("confidence"), level = 0.50,type="response") 
polygon(c(rev(prep.grd), prep.grd), 
		c(rev(prd[,2]), prd[,3]), 
        col=rgb(0, 0, 0.5,0.25),
		border = NA)
	
x1 = c(950-90, 950+90)
y1 = c(500,500)
y2 = c(-500,-500)	

abline(v = 950, col = "blue", lwd = 4, lty = 2)

polygon(c(x1, rev(x1)), c(y1, rev(y2)),
      col=rgb(0, 0, 0.5,0.25), 
border = NA)

# overlay trendy	
d.temporal.trendy = d.temporal.trendy[order(d.temporal.trendy$prep),]
mylm = 	lm(mean~prep, data = d.temporal.trendy)
abline(mylm,col=rgb(178/255, 178/255, 0,1), lwd = 3)

prd<-predict(mylm,newdata=data.frame(prep = prep.grd),interval = c("confidence"), level = 0.50,type="response") 
polygon(c(rev(prep.grd), prep.grd), 
		c(rev(prd[,2]), prd[,3]), 
        col=rgb(178/255, 178/255, 0,0.25),
		border = NA)
		
legend("bottomleft", 
       # inset=0.05, 
	   # legend = c("ENF","DBF","MF","SHB","GRA","CRO"),
#	   legend = c("EC", "RS","TRENDY"),
       	   legend = c("Constrained EC obs", "Constrained Global Obs","TRENDY"),
	   horiz=F,
	   lwd = 4,
	   col = c(rgb(0, 1, 0,1), rgb(0, 0, 1,1), rgb(178/255, 178/255, 0,1)),
	   text.col = c(rgb(0, 1, 0,1), rgb(0, 0, 1,1), rgb(178/255, 178/255, 0,1)),
	   cex = 1.5,
	   box.col = "transparent",
       bg = "transparent")	

text(x = 100, y = 95, "a)",cex = 2)

# abline(v = c(650, 750, 700 + 100), lty = c(2,1,2))
# abline(v = c(650, 700 + 100), lty = c(2,2))

# polygon(c(x1, rev(x1)), c(y1, rev(y2)),
      # col = rgb(0.5,0.5,0.5,0.3), 
# border = NA)
	   
## plot spatial sensitivity
		
# par(mar=c(0,5,0,0))
plot(mean~prep, data = d.spatial.ec, type="l", 
								  ylim = c(-100,100),
                                  xlim = c(100, 1300),								  
								  cex = 2, lwd = 4, 
								  col = "green",
								  xlab = "Anuual Mean Precipitation (MAP: mm)", 
								  # ylab = "Spatial sensitivity", 
								  # ylab = expression(paste(beta ["spatial"])),
								  ylab = "",
								  cex.axis = 1.5, cex.lab = 1.6,yaxt='n')

# add axis and labels
axis(side = 2, las = 1, tck = -.02, at = c(-50,0,50), labels = c("-50","0","50"),cex.axis = 1.5)


polygon(c(rev(d.spatial.ec$prep), d.spatial.ec$prep), 
		c(rev(d.spatial.ec$mean-d.spatial.ec$sd), d.spatial.ec$mean+d.spatial.ec$sd), 
        col=rgb(0, 0.5, 0,0.25),
		border = NA)						  
								  
abline(h = 0, lty = 2, lwd = 2)		

x1 = c(750-75,750+60)
y1 = c(500,500)
y2 = c(-500,-500)
abline(v = 750, col = "green", lwd = 4,lty = 2)
# abline(v = c(650, 700 + 100), lty = c(2,2))

polygon(c(x1, rev(x1)), c(y1, rev(y2)),
      col=rgb(0, 0.5, 0,0.25), 
border = NA)

# overlay remote sensing
par(new=TRUE)
plot(mean~prep, data = d.spatial.rs, type="l", 
								  ylim = c(-100,100),
                                  xlim = c(100, 1300),								  
								  col="blue",lwd = 4,
								  bty='n',pch='',ylab='',xlab='',yaxt='n',xaxt='n', ann=FALSE)

polygon(c(rev(d.spatial.rs$prep), d.spatial.rs$prep), 
		c(rev(d.spatial.rs$mean-d.spatial.rs$sd), d.spatial.rs$mean+d.spatial.rs$sd), 
 col=rgb(0, 0, 0.5,0.25),
 border = NA)
	
x1 = c(830-70, 830+70)
y1 = c(500,500)
y2 = c(-500,-500)	

polygon(c(x1, rev(x1)), c(y1, rev(y2)),
      col=rgb(0, 0, 0.5,0.25), 
border = NA)

abline(v = 830, col = "blue", lwd = 4, lty = 2)	

# overlay trendy	
par(new=TRUE)
plot(mean~prep, data = d.spatial.trendy, type="l", 
								  ylim = c(-100,100),
                                  xlim = c(100, 1300),								  
								  col=rgb(178/255, 178/255, 0,1),
								  lwd = 4,
								  bty='n',pch='',ylab='',xlab='',yaxt='n',xaxt='n', ann=FALSE)

		polygon(c(rev(d.spatial.trendy$prep), d.spatial.trendy$prep), 
		c(rev(d.spatial.trendy$mean-d.spatial.trendy$sd), d.spatial.trendy$mean+d.spatial.trendy$sd), 
        col=rgb(178/255, 178/255, 0,0.25),
        border = NA)

text(x = 100, y = 95, "b)",cex = 2)

# x1 = c(650,800)
# y1 = c(500,500)
# y2 = c(-500,-500)
# abline(v = c(650, 750, 700 + 100), lty = c(2,1,2))
# abline(v = c(650, 700 + 100), lty = c(2,2))


# polygon(c(x1, rev(x1)), c(y1, rev(y2)),
      # col = rgb(0.5,0.5,0.5,0.3), 
# border = NA)

	   
mtext(side = 1, line = 3, 
      "Mean Anuual Precipitation (MAP: mm)", 
     # "mm",
      outer = TRUE, cex = 1.6, col = rgb(0, 0, 0,1))
				
# mtext(side = 2, line = 3, 
#      expression("" ~ g ~ C ~ m^{-2} ~ yr ^{-1}~ "per 100mm"), 
#      outer = TRUE, cex = 1.6, col = rgb(0, 0, 0,1))

# add spatial sensitivity difference
mtext(side = 2, line = 3, 
      expression("" ~ Delta * delta ^{s} ~ ""), 
      outer = TRUE, cex = 1.6, adj = 0.23, col = rgb(0, 0, 0,1))

# add temproal sensitivity difference
mtext(side = 2, line = 3, 
      expression("" ~ Delta * delta ^{t} ~ ""), 
      outer = TRUE, cex = 1.6, adj = 0.78, col = rgb(0, 0, 0,1))

dev.off()		
		
		

####new figure 2-- for temperature
#plot together
d.temporal.ec = data.frame(read.csv("F:/zhihua/dataset/results2/delt.temporal.ec.temp.csv"))[,-1]
d.spatial.ec = data.frame(read.csv("F:/zhihua/dataset/results2/delt.spatial.ec.temp.csv"))[,-1]

d.temporal.rs = data.frame(read.csv("F:/zhihua/dataset/results2/delt.temporal.rs.temp.csv"))[,-1]
d.spatial.rs = data.frame(read.csv("F:/zhihua/dataset/results2/delt.spatial.rs.temp.csv"))[,-1]

d.temporal.trendy = data.frame(read.csv("F:/zhihua/dataset/results2/delt.temporal.trendy.temp.csv"))[,-1]
d.spatial.trendy = data.frame(read.csv("F:/zhihua/dataset/results2/delt.spatial.trendy.temp.csv"))[,-1]

prep.grd = seq(1,25, length.out = 25)

png("F:/zhihua/dataset/results2/fig2-3-temperature.png",height = 2500, width = 2000, res = 300, units = "px")

par(mfrow=c(2,1),mar=c(0,0,0,0),oma=c(5,5,0,0))

 # plot temporal sensitivity
# http://htmlcolorcodes.com/
# par(mar=c(0,0,0,0))

plot(1, type="n", 
								  ylim = c(-100,100),
                                  xlim = c(1, 25),	
                                  # xlim = c(100, 1300),								  
								  cex = 2, lwd = 4, 
								  col = "green",
								  xlab = "Anuual Mean Precipitation (MAP: mm)", 
								  # ylab = "Spatial sensitivity", 
								  # ylab = expression(paste(beta ["temporal"])), 
								  cex.axis = 1.5, 
								  cex.lab = 1.6,
								  xaxt='n', ann=FALSE,yaxt='n')

# add axis and labels
axis(side = 2, las = 1, tck = -.02, at = c(-50,0,50), labels = c("-50","0","50"), cex.axis = 1.5)

abline(h = 0, lty = 2, lwd = 2)		

points(x = d.temporal.ec$prep, y = d.temporal.ec$mean, pch = 0, cex = 2,col="green")
points(x = d.temporal.rs$prep, y = d.temporal.rs$mean, pch = 1, cex = 2,col="blue")
points(x = d.temporal.trendy$prep, y = d.temporal.trendy$mean, pch = 2, cex = 2,col=rgb(178/255, 178/255, 0,1))

text(x = 1.5, y = 95, "a)",cex = 2)

## plot spatial sensitivity	
plot(mean~prep, data = d.spatial.ec, type="l", 
								  ylim = c(-100,100),
                                  xlim = c(1, 25),								  
								  cex = 2, lwd = 4, 
								  col = "green",
								  xlab = "Anuual Mean Temperature (degree)", 
								  # ylab = "Spatial sensitivity", 
								  # ylab = expression(paste(beta ["spatial"])),
								  ylab = "",
								  cex.axis = 1.5, cex.lab = 1.6,yaxt='n')

# add axis and labels
axis(side = 2, las = 1, tck = -.02, at = c(-50,0,50), labels = c("-50","0","50"), cex.axis = 1.5)

polygon(c(rev(d.spatial.ec$prep), d.spatial.ec$prep), 
		c(rev(d.spatial.ec$mean-d.spatial.ec$sd), d.spatial.ec$mean+d.spatial.ec$sd), 
        col=rgb(0, 0.5, 0,0.25),
		border = NA)						  
								  
abline(h = 0, lty = 2, lwd = 2)		

x1 = c(750-75,750+60)
y1 = c(500,500)
y2 = c(-500,-500)
abline(v = 750, col = "green", lwd = 4,lty = 2)
# abline(v = c(650, 700 + 100), lty = c(2,2))

polygon(c(x1, rev(x1)), c(y1, rev(y2)),
      col=rgb(0, 0.5, 0,0.25), 
border = NA)

# overlay remote sensing
par(new=TRUE)
plot(mean~prep, data = d.spatial.rs, type="l", 
								  ylim = c(-100,100),
                                  xlim = c(1, 25),								  
								  col="blue",lwd = 4,
								  bty='n',pch='',ylab='',xlab='',yaxt='n',xaxt='n', ann=FALSE)

polygon(c(rev(d.spatial.rs$prep), d.spatial.rs$prep), 
		c(rev(d.spatial.rs$mean-d.spatial.rs$sd), d.spatial.rs$mean+d.spatial.rs$sd), 
 col=rgb(0, 0, 0.5,0.25),
 border = NA)
	
x1 = c(830-70, 830+70)
y1 = c(500,500)
y2 = c(-500,-500)	

polygon(c(x1, rev(x1)), c(y1, rev(y2)),
      col=rgb(0, 0, 0.5,0.25), 
border = NA)

abline(v = 830, col = "blue", lwd = 4, lty = 2)	

# overlay trendy	
par(new=TRUE)
plot(mean~prep, data = d.spatial.trendy, type="l", 
								  ylim = c(-100,100),
                                  xlim = c(1, 25),								  
								  col=rgb(178/255, 178/255, 0,1),
								  lwd = 4,
								  bty='n',pch='',ylab='',xlab='',yaxt='n',xaxt='n', ann=FALSE)

		polygon(c(rev(d.spatial.trendy$prep), d.spatial.trendy$prep), 
		c(rev(d.spatial.trendy$mean-d.spatial.trendy$sd), d.spatial.trendy$mean+d.spatial.trendy$sd), 
        col = rgb(178/255, 178/255, 0,0.25),
        border = NA)

text(x = 1.5, y = 195, "b)",cex = 2)

legend("bottomleft", 
       	   legend = c("Constrained EC obs", "Constrained Global Obs","TRENDY"),
	   horiz=F,
	   # pch = 0:2,
	   lwd = 4,
	   col = c(rgb(0, 1, 0,1), rgb(0, 0, 1,1), rgb(178/255, 178/255, 0,1)),
	   text.col = c(rgb(0, 1, 0,1), rgb(0, 0, 1,1), rgb(178/255, 178/255, 0,1)),
	   cex = 1.5,
	   box.col = "transparent",
       bg = "transparent")	

# polygon(c(x1, rev(x1)), c(y1, rev(y2)),
      # col = rgb(0.5,0.5,0.5,0.3), 
# border = NA)

	   
mtext(side = 1, line = 3, "Anuual Mean Temperature (degree)", 
      outer = TRUE, cex = 1.6, col = rgb(0, 0, 0,1))
				
# mtext(side = 2, line = 3, 
      # expression("" ~ g ~ C ~ m^{-2} ~ yr ^{-1}~ " per degree"), 
 #    expression("" ~ g ~ C ~ m^{-2} ~ yr ^{-1}~ ""), 
 #     outer = TRUE, cex = 1.6, col = rgb(0, 0, 0,1))

# add spatial sensitivity difference
mtext(side = 2, line = 3, 
      expression("" ~ Delta * gamma ^{s} ~ ""), 
      outer = TRUE, cex = 1.6, adj = 0.23, col = rgb(0, 0, 0,1))

# add temproal sensitivity difference
mtext(side = 2, line = 3, 
      expression("" ~ Delta * gamma ^{t} ~ ""), 
      outer = TRUE, cex = 1.6, adj = 0.78, col = rgb(0, 0, 0,1))


dev.off()		
		

