required <- c('cowplot','emmeans','ggplot2','influence.ME','lme4','lmerTest','magick','pROC','performance','survival','pbkrtest')
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop('Install required packages: ',paste(missing,collapse=', '))
print(data.frame(Package=required,Version=vapply(required,function(p)as.character(packageVersion(p)),character(1))),row.names=FALSE)
cat('Dependency preflight passed. Kenward-Roger inference requires pbkrtest.\n')
