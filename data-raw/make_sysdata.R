# Assemble R/sysdata.rda from the calibration objects. Run from data-raw/
# after the scripts that produce them:
#   dfgls_rs    Rscript dfgls_response_surfaces.R dfgls_rs.rds   (D11)
#   kpss_limit  Rscript kpss_limit.R                              (D14)
#   phi_rs, rho_rs  Rscript new_response_surfaces.R                (D16)
dfgls_rs <- readRDS("dfgls_rs.rds")
kpss_limit <- readRDS("kpss_limit.rds")
phi_rs <- readRDS("phi_rs.rds")
rho_rs <- readRDS("rho_rs.rds")
objs <- c("dfgls_rs", "kpss_limit", "phi_rs", "rho_rs")
save(list = objs, file = "../R/sysdata.rda", compress = "xz")
cat("R/sysdata.rda:", objs, "\n")
