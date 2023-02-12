temp <-
  gen_fields() %>%
assign_trial_design(., design = c("Latin Square Fixed 5", "Latin Square Fixed 6")) %>%
  gen_field_parameters(., nsim = 2) %>%
#--- generate trial design ---#
dplyr::mutate(trial_design = list(
  gen_trial_design(field_sf, field_pars, design_name, num_treatments)
))

temp$trial_design
