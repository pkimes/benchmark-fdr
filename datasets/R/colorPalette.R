# Custom color palette to be consistent across plots and case studies
# Specific color choices adapted from Alyssa Frazee's RSkittleBrewer 
# package (https://github.com/alyssafrazee/RSkittleBrewer)
# Also including ideas from discussion in meeting on 4/6/2018

# Bonferroni - dashed black
# BH - solid black
# qvalue - gray45
# ashq - red3
# ashs - hotpink
# ihw - green3
# bl - darkorange1
# locfdr - dodgerblue3
# scott empirical - purple4
# scott theoretical - mediumpurple3
# AdaPT GLM - darkgoldenrod
# AdaPT GAM - darkgoldenrod3

candycols <- data.frame(Method = c("bonf", "bh", "qvalue", "ashs", "ashq", 
                                   "ihw", "bl", "lfdr", 
                                   "fdrreg-t", "fdrreg-e",
                                   "adapt-glm", "adapt-gam"),
                        col = c("black", "black", "gray45", "hotpink", "red3",
                                "green3", "darkorange1", "dodgerblue3", 
                                "purple4", "mediumpurple3",
                                "darkgoldenrod", "darkgoldenrod3"),
                        lty = c("dashed", rep("solid", 11)),
                        stringsAsFactors = FALSE)
