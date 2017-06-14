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

color1 = c("grey20", "grey90")
ylab=expression("" ~ g ~ C ~ m^{-2} ~ yr ^{-1}~ " Per 100 mm")

p1 <- ggplot(data=D1, aes(x=region, y=mean, fill=flux)) + 
    geom_bar(colour="black", stat="identity",
             position=position_dodge(),
             size=.3) + 
    facet_grid(type~Method)+			 
	ylab(ylab) + 	 
    # ylab(expression(paste(beta ["temporal"]))) + 
	# xlab("Region by MAP") + 
	# ylab("") + 
	xlab("Region By Map (mm)") + 
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
    theme(strip.text.y = element_text(size=18)) +
    scale_fill_manual(values=color1, 
                    name="",
                    breaks=levels(D1$flux),
                    labels=levels(D1$flux)) 

print(p1)

ggsave("F:/zhihua/dataset/results2/sensitivity.png", width = 8, height = 6, units = "in")

