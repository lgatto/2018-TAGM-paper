numMarkers <- sum(fData(andyMCMChl)$markers != "unknown")
numProteins <- nrow(andyMCMChl)
svmFDR <- sum(fData(andyMCMChl)$final.assignment!= "unknown")
numMCMC <- sum((as.numeric(fData(andyMCMChl)$predicted.probability) * fData(andyMCMChl)$outlier.allocation) > 0.99)
numShannon <- sum(fData(andyMCMChl)$meanshannon < 0.1)

variables <- c("Input","Input","SVM", "SVM", "SVM", "TAGM", "TAGM", "TAGM", "TAGM")
values <- c(numMarkers, numProteins - numMarkers, numMarkers, svmFDR-numMarkers, numProteins - svmFDR,
            numMarkers, numMCMC  - numMarkers, numShannon - numMCMC, numProteins - numShannon)

values2 <- factor(c("Markers", "Unknown", "Markers", "<5% FDR", "Unknown",
                    "Markers",">0.99 Posterior probability", "<0.1 Shannon Entropy", "Unknown"),
                  levels = rev(c("Markers","<5% FDR",">0.99 Posterior probability", "<0.1 Shannon Entropy", "Unknown")))
df <- data.frame(variables, values, values2)

gg <- ggplot(df, aes(variables, values))
gg <- gg + geom_bar(stat = "identity", aes(fill= values2), alpha = 0.5) + scale_fill_manual(values = rev(c("blue", "limegreen", " limegreen","coral","light blue")))
gg <- gg + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
gg <- gg + annotate("text", x = c(1,2,3), y = 500, label = "Markers", size = 10) + theme(legend.position="none")
gg <- gg + annotate("text", x = c(1,2,3), y = 4500, label = "Unknown", size = 10)
gg <- gg + annotate("text", x = 2, y = 2000, label = "<5% FDR", size = 8)
gg <- gg + annotate("text", x = 3, y = 1800, label = "Posterior Probability", size = 8)
gg <- gg + annotate("text", x = 3, y = 2200, label = ">0.99" , size = 8)
gg <- gg + annotate("text", x = 3, y = 3800, label = "<0.1" , size = 8)
gg <- gg + annotate("text", x = 3, y = 3400, label = "Shannon Entropy", size = 8)
gg <- gg + ggtitle("Effect of Methodology on Protein Assignment") +
  xlab("Method") + ylab("Number of Proteins") +  theme(plot.title = element_text(hjust = 0.5))
gg
