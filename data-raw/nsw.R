
## The original data, nsw.dta, is the NSW dataset analyzed by Smith and Todd (2005).
# The file was downloaded from
## https://dataverse.harvard.edu/file.xhtml?persistentId=doi:10.7910/DVN/23407/DYEWLO&version=1.0.


# Data in "wide format"
nsw <- haven::read_dta("nsw.dta")
nsw <- as.data.frame(nsw)

#Create "experimental" treatment dummy: 1 if in experimental sample, 0 if in non-experimental
nsw$experimental <- ifelse(is.na(nsw$treated), 0 , 1)
usethis::use_data(nsw, overwrite = TRUE)


# Data in "long format"
nsw_long <- nsw
nsw_long$id <- 1:dim(nsw)[1]
nsw_long <- panelr::long_panel(nsw_long, id = "id",
                         wave = "year",
                         prefix = "",
                         periods = c(75,78),
                         label_location = "end")
nsw_long$year = ifelse(nsw_long$year==75, 1975, 1978)
nsw_long <- as.data.frame(nsw_long)
usethis::use_data(nsw_long, overwrite = TRUE)


