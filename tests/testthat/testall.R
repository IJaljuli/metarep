data(CD007077_CMP001)
CD007077 <- CD007077_CMP001
# data.frame(n.e = rep( 25 , 7 ) ,
#                             n.c = rep( 25 , 7 ) ,
#                             mean.e = rep( 0, times = 7 ),
#                             mean.c = rnorm(n = 7,mean = 1),
#                             sd.e = round(sqrt(rchisq(n = 7,df = 24)/24),digits=2),
#                             sd.c = round(sqrt(rchisq(n = 7,df = 24)/24),digits=2) )

m1 <- meta::metabin( event.e = N_EVENTS1, n.e = N_TOTAL1, 
                     event.c = N_EVENTS2, n.c = N_TOTAL2,
                     studlab = STUDY, comb.fixed = T , comb.random = F,
                     method = 'MH', sm = CD007077$SM[1],
                     data = CD007077)
m1.ra <- metarepl(x = m1 , u = 2 , common.effect = F ,t = 0.05 ,report.u.max = F)
m1.ra.bounds <- metarepl(x = m1 , u = 2 , common.effect = F ,t = 0.05 ,report.u.max = T)
find_umax(x = m1 , common.effect = F,alternative = 'less',t = 0.05,confidence = 0.975)
find_umax(x = m1 , common.effect = F,alternative = 'two-sided',t = 0.05,confidence = 0.95)
