# Load data
filePath <- "/Users/woodie/Documents/Courses/ISyE 6416 Computational Statistics (Spring 2018)/HW/ISYE-6416/hw1/w_logret_3automanu.csv"
rawdata  <- read.csv(file=filePath, header=FALSE, sep=",")
colnames(rawdata) <- c("Toyota Motor Corp.", "Ford Motor Corp.", "GM.")
summary(rawdata)

# Plot the 3d surface for response (GM.)
plot_ly(x=rawdata$`Toyota Motor Corp.`,
        y=rawdata$`Ford Motor Corp.`,
        z=rawdata$`GM.`, type="mesh3d")

lmResult1 <- lm(`GM.` ~ `Ford Motor Corp.` + `Toyota Motor Corp.`, # regression formula
                data=rawdata)                                      # data set
summary(lmResult1)

lmResult2 <- lm(`GM.` ~ `Ford Motor Corp.`, # regression formula
                data=rawdata)               # data set
summary(lmResult2)

lmResult3 <- lm(`GM.` ~ `Toyota Motor Corp.`, # regression formula
                data=rawdata)                 # data set
summary(lmResult3)

AIC(lmResult1)
BIC(lmResult1)