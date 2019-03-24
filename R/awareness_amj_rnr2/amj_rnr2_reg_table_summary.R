library(xml2)
library(rvest)


path <- 'C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/compnet2/amj_rnr2_results/firmctrl'

R <- 2000
nPeriod <- 11
mx <- 'm4'

items <- c(
  'clarabridge',
  'confirmit',
  'getfeedback',
  'medallia',
  'snap-surveys-ltd',
  'checkmarket',
  'cloudcherry',
  'customergauge',
  'cx-index',
  'empathica',
  'feedback-lite',
  'first-mile-geo',
  'inqwise',
  'leaderamp',
  'myfeelback',
  'promoter-io',
  'satmetrix',
  'super-simple-survey',
  'survata',
  'surveygizmo',
  'surveyrock',
  'typeform',
  'userate',
  'verint',
  'voice-polls'
)

df <- data.frame()
for (name in items)
{
  cat(sprintf(' %s\n',name))
  htmlfile <- sprintf('%s_tergm_results_pd%s_R%s_%s.html', name, nPeriod, R, mx)
  file_path <- file.path(path,htmlfile)
  
  if (!file.exists(file_path))
    next
  
  dat <- read_html(file_path)
  node <- html_node(dat, 'table')
  reg <- html_table(node, fill=T)
  
  j.name <- which(names(reg)=="")
  j.mod <- which( ! names(reg) %in% c(NA,'NA',''))
  
  i.ci <- which(reg[,j.name]=="")
  i.eff <- which( ! reg[,j.name] %in% c("","* 0 outside the confidence interval"))
  
  if (length(names(df))==0)  {
    df <- data.frame(.=reg[i.eff, j.name])    
  }
  
  tmpdf <- data.frame(a=reg[i.eff, j.mod], b=c(reg[i.ci, j.mod],""))
  names(tmpdf) <- c(name, '.')

  df <- cbind(df, tmpdf)
}

