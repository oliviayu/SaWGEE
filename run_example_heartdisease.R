library(magrittr)
library(dplyr)

source("sawgee_code.r")

df=read.csv("heart_2022_no_nans.csv")
xvars = setdiff(names(df), c("HadHeartAttack"))
df$yvar = as.numeric(df$HadHeartAttack=='Yes')

df2 = data.frame(idvar=1:length(unique(df$State)), State=unique(df$State)) %>%
  left_join(df, ., by='State') %>% arrange(idvar)
df2$Sex = as.numeric(df2$Sex=='Male')
df2$HadAngina = as.numeric(df2$HadAngina=='Yes')
df2$HadStroke = as.numeric(df2$HadStroke=='Yes')
df2$HadDiabetes = as.numeric(df2$HadDiabetes=='Yes')
df2$BlindOrVisionDifficulty = as.numeric(df2$BlindOrVisionDifficulty=='Yes')
df2$DifficultyWalking = as.numeric(df2$DifficultyWalking=='Yes')
df2$ChestScan = as.numeric(df2$ChestScan=='Yes')
df2$AlcoholDrinkers = as.numeric(df2$AlcoholDrinkers=='Yes')
df2$FluVaxLast12 = as.numeric(df2$FluVaxLast12=='Yes')
df2$GeneralHealth=as.factor(df2$GeneralHealth)
df2$lct=as.numeric(df2$LastCheckupTime == "Within past year (anytime less than 12 months ago)")
df2$RemovedTeeth = as.factor(df2$RemovedTeeth)
df2$smoke= as.numeric(df2$SmokerStatus %in% c("Current smoker - now smokes every day", "Current smoker - now smokes some days"))
df2$SleepHours = scale(df2$SleepHours)
df2$HeightInMeters = scale(df2$HeightInMeters)
df2$age = NA
df2$age[df2$AgeCategory %in% c("Age 18 to 24", "Age 25 to 29", "Age 30 to 34", 
                               "Age 35 to 39")] = "Age 18 to 39"
df2$age[df2$AgeCategory %in% c("Age 40 to 44", "Age 45 to 49", "Age 50 to 54")] = 'Age 40 to 54'
df2$age[df2$AgeCategory %in% c("Age 55 to 59", "Age 60 to 64", "Age 65 to 69")] = 'Age 55 to 69'
df2$age[df2$AgeCategory %in% c("Age 70 to 74", "Age 75 to 79", "Age 80 or older")] = 'Age 70 or older'

sxvars = c("Sex","GeneralHealth","lct","SleepHours",
           "RemovedTeeth","HadAngina", "HadStroke",
           "HadDiabetes", "BlindOrVisionDifficulty",
           "DifficultyWalking", "smoke","ChestScan", 
           "age","AlcoholDrinkers","FluVaxLast12","HeightInMeters")

fm = as.formula(paste("yvar~",paste(sxvars, collapse = '+')))
nvmd1 = glm(fm, family = 'binomial', data=df2)
summary(nvmd1)

mu = predict.glm(nvmd1, df2, type="response")
ini.beta = as.numeric(coef(nvmd1))
rho_ini  = 0

mean(df$yvar)*2
frac=0.11
M = round(frac*nrow(df2)/2)*2

T_n = 150
parallel=F; mc.cores=1

### Hybrid with Lopt criteria
OS.result_hybrid = subWGEE(mu, df2, M, T_n, ini.beta, fm, strategy = 'os',
                           hybrid = T,
                           parallel = parallel, mc.cores = mc.cores, 
                           method = 1,  rho_ini=rho_ini)
summary_subWGEE(OS.result_hybrid$output[,-1], as.numeric(ini.beta))

### Lopt criteria
OS.result = subWGEE(mu, df2, M, T_n, ini.beta, fm, strategy = 'os',
                    hybrid = F,
                    parallel = parallel, mc.cores = mc.cores, 
                    method = 1,  rho_ini=rho_ini)
summary_subWGEE(OS.result$output[,-1], as.numeric(ini.beta))

## random sampling
RS.result = subWGEE(mu, df2, M, T_n, ini.beta, fm, strategy = 'random',
                    parallel = parallel, mc.cores = mc.cores)
summary_subWGEE(RS.result$output[,-1], as.numeric(ini.beta))
