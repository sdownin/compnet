library(igraph)
library(network)
library(intergraph)
library(stringr)
library(plyr)

setwd("C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\compnet2\\firm_nets_rnr")

skip <- c('surveygizmo_d3.rds')

files <- dir(pattern = "d3\\.rds$")
files <- files[!files %in% skip]

cdf <- data.frame()

for (file in files) {
  
  focal <- str_split(file, "[_]")[[1]][1]
  cat(sprintf('loading firm %s file %s\n', focal, file))
  
  nets <- readRDS(file)
  gs <- lapply(nets, asIgraph)
  g <- gs$`2017`
  
  df <- data.frame(
    firm = V(g)$vertex.names,
    domain = V(g)$domain,
    focal = focal,
    d = as.character(c(igraph::distances(g, which(V(g)$vertex.names==focal), V(g)))),
    stringsAsFactors = F
  )
  
  df$d[ which(!df$d %in% c('0','1','2','3','4')) ] <- '5+'
  df <- df[order(df$d, decreasing = F), ]
  
  cdf <- rbind(cdf, df)
  
}

all.focal.firms <- unique(cdf$focal)

ucdf <- unique(cdf[,1:2])

cnt <- plyr::count(cdf$firm)
cnt <- cnt[order(cnt$freq, decreasing = T),]

head(cnt, 50)

focal.firms <- c(
  'qualtrics',
  'clarabridge','confirmit','medallia','snap-surveys-ltd','getfeedback' ## empathica renamed InMoment 
)

cdf2 <- cdf[which(cdf$focal %in% focal.firms),]

cnt2 <- plyr::count(cdf2$firm)
cnt2 <- cnt2[order(cnt2$freq, decreasing = T),]

cdf2r <- ddply(cdf2, .(firm),  summarize,
               domain=unique(domain),
               nnets=length(focal),
               focal=paste(focal, collapse="|"))
cdf2r <- cdf2r[order(cdf2r$nnets, decreasing = T), ]

out.dir <- "C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\compnet2"
write.csv(cdf2r, file = file.path(out.dir, "six_focal_firm_unique_indirect_comp_net_firms.csv"), row.names = F)

## add employee count
# g.full.file <- "C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\compnet2\\g_full.graphml"
# g.full <- read.graph(g.full.file, format='graphml')
# vdf <- igraph::as_data_frame(g.full, what='vertices')
codf.file <- "C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\crunchbase\\crunchbase_export_20161024\\organizations.csv"
co.df <- read.csv(codf.file, header = T, stringsAsFactors = F)
co.df$geo <- apply(co.df[ ,c('country_code','state_code','city')],1,function(x)paste(x,collapse=", "))
cdf2e <- merge(cdf2r, co.df[,c('company_name_unique','employee_count','status','geo')], by.x='firm',by.y='company_name_unique',all.x=T,all.y=F)

out.dir <- "C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\compnet2"
write.csv(cdf2e, file = file.path(out.dir, "six_focal_firm_unique_indirect_comp_net_firms_WEMPLOYEES.csv"), row.names = F)

ORIG_6_NET_FIRMS <- cdf2e

# library(curl)
library(jsonlite)
library(rvest)
library(stringr)
library(urltools)
library(curl)

##
#  fetch current/active url for old domain
##
fetchUrlSafe <- function(url, domain.only=T) {
  req <- tryCatch(
    curl_fetch_memory(url),
    error = function(e) e
  )
  if(inherits(req,  "error"))
    return("ERROR")
  if (domain.only)
    return(gsub("www\\.","",domain(req$url)))
  return(req$url)
}

## fetch updated urls
cdf2ed <- cdf2e
cdf2ed$new_domain <- NA
for (i in 1:nrow(cdf2ed)) {
  if (i %% 2 == 0) 
    cat(sprintf("i = %s (%5.1f%s)  %s\n",i,100*i/nrow(cdf2ed),'%',cdf2ed$domain[i]))
  cdf2ed$new_domain[i] <- fetchUrlSafe(cdf2ed$domain[i])
}

out.dir <- "C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\compnet2"
write.csv(cdf2ed, file = file.path(out.dir, "six_focal_firm_unique_indirect_comp_net_firms_WEMPLOYEES_NEWDOMAIN.csv"), row.names = F)



##------- OTHER FOCAL FIRMS AFTER FIRST 6  -----------------------------------

# focal.firms <- c(
#   'qualtrics',
#   'clarabridge','confirmit','medallia','snap-surveys-ltd','getfeedback' ## empathica renamed InMoment 
# )

focal.firms <- c(
  ## D1
  'checkmarket','cloudcherry','cx-index','feedback-lite','first-mile-geo',
  'leaderamp','myfeelback','promoter-io','super-simple-survey',
  'survata','surveyrock','typeform','userate',
  ## D2
  'customergauge','inqwise','verint','voice-polls',
  ## D3
  'empathica','satmetrix'
)

cdf2 <- cdf[which(cdf$focal %in% focal.firms),]

cnt2 <- plyr::count(cdf2$firm)
cnt2 <- cnt2[order(cnt2$freq, decreasing = T),]

cdf2r <- ddply(cdf2, .(firm),  summarize,
               domain=unique(domain),
               nnets=length(focal),
               focal=paste(focal, collapse="|"))
cdf2r <- cdf2r[order(cdf2r$nnets, decreasing = T), ]

out.dir <- "C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\compnet2"
write.csv(cdf2r, file = file.path(out.dir, "19_focal_firm_unique_indirect_comp_net_firms.csv"), row.names = F)

## add employee count
# g.full.file <- "C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\compnet2\\g_full.graphml"
# g.full <- read.graph(g.full.file, format='graphml')
# vdf <- igraph::as_data_frame(g.full, what='vertices')
codf.file <- "C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\crunchbase\\crunchbase_export_20161024\\organizations.csv"
co.df <- read.csv(codf.file, header = T, stringsAsFactors = F)
co.df$geo <- apply(co.df[ ,c('country_code','state_code','city')],1,function(x)paste(x,collapse=", "))
cdf2e <- merge(cdf2r, co.df[,c('company_name_unique','employee_count','status','geo')], by.x='firm',by.y='company_name_unique',all.x=T,all.y=F)

out.dir <- "C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\compnet2"
write.csv(cdf2e, file = file.path(out.dir, "19_focal_firm_unique_indirect_comp_net_firms_WEMPLOYEES.csv"), row.names = F)

## limit to firms not in original 6 focal firm networks
cdf2eno <- cdf2e[ ! cdf2e$firm %in% ORIG_6_NET_FIRMS$firm, ]

## fetch updated urls
cdf2ed <- cdf2eno
cdf2ed$new_domain <- NA
for (i in 1:nrow(cdf2ed)) {
  if (i %% 10 == 0) 
    cat(sprintf("i = %s (%5.1f%s)  %s\n",i,100*i/nrow(cdf2ed),'%',cdf2ed$domain[i]))
  cdf2ed$new_domain[i] <- fetchUrlSafe(cdf2ed$domain[i])
}

out.dir <- "C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\compnet2"
write.csv(cdf2ed, file = file.path(out.dir, "19_focal_firm_unique_indirect_comp_net_firms_WEMPLOYEES_NEWDOMAIN.csv"), row.names = F)



## add closed on dates
.closed <- co.df[,c('company_name_unique','closed_on')]
cdf2edc <- merge(cdf2ed, co.df[,c('company_name_unique','closed_on')], by.x='firm',by.y='company_name_unique',all.x=T,all.y=F)

out.dir <- "C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\compnet2"
write.csv(cdf2edc, file = file.path(out.dir, "19_focal_firm_CLOSED_ON.csv"), row.names = F, na = "-")



##------ QUALTRICS COMPETITORS BY DISTANCE -----------------------------------
nets <- readRDS("qualtrics_d3.rds")
gs <- lapply(nets, asIgraph)
g <- gs$`2017`

df <- data.frame(
  firm = V(g)$vertex.names,
  domain = V(g)$domain,
  d = as.character(c(igraph::distances(g, which(V(g)$vertex.names=='qualtrics'), V(g)))),
    stringsAsFactors = F
)

df$d[ which(!df$d %in% c('0','1','2','3','4')) ] <- '5+'
df <- df[order(df$d, decreasing = F), ]

out.dir <- "C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\compnet2"
write.csv(df, file = file.path(out.dir, "qualtrics_d3_2016_comp_dist_table.csv"), row.names = F)
