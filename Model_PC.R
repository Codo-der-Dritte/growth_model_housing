bs_pc <- sfcr_matrix(
  columns = c("Households C", "Households W", "Firms", "Government", "Central bank", "sum"),
  codes = c("hc", "hw",  "f", "g", "cb", "s"),
  r1 = c("Money", hc = "+Hhc", hw = "+Hhw", cb = "-Hs"),
  r2 = c("Bills", hc = "+Bh", g = "-Bs", cb = "+Bcb"),
  r3 = c("Balance", hc = "-Vc", hw = "-Vw", g = "+V")
)

sfcr_matrix_display(bs_pc, "bs")

tfm_pc <- sfcr_matrix(
  columns = c("Households C", "Households W",  "Firms", "Government", "CB current", "CB capital"),
  codes = c("hc", "hw",  "f", "g", "cbc", "cbk"),
  c("Consumption", hc = "-Cc", hw = "-Cw", f = "+C"),
  c("Govt. Expenditures", f = "+G", g = "-G"),
  c("Income", hc = "+Yc", hw = "+Yw", f = "-Y"),
  c("Int. payments", hc = "+r[-1] * Bh[-1]", g = "-r[-1] * Bs[-1]", cbc = "+r[-1] * Bcb[-1]"),
  c("CB profits", g = "+r[-1] * Bcb[-1]", cbc = "-r[-1] * Bcb[-1]"),
  c("Taxes", hc = "-TXc", hw = "-TXw",  g = "+TX"),
  c("Ch. Money", hc = "-(Hhc - Hhc[-1])", hw = "-(Hhw - Hhw[-1])", cbk = "+(Hs - Hs[-1])"),
  c("Ch. Bills", hc = "-(Bh - Bh[-1])", g = "+(Bs - Bs[-1])", cbk = "-(Bcb - Bcb[-1])")
)

sfcr_matrix_display(tfm_pc)

pc_eqs <- sfcr_set(
  Yc ~ Cc + G * zeta,
  Yw ~ Cw + G * (1-zeta),
  Y ~ Yc + Yw,
  YDc ~ Yc - TXc + r[-1] * Bh[-1],
  YDw ~ Yw - TXw,
  YD ~ YDc + YDc,
  TXc ~ theta * (Yc + r[-1] * Bh[-1]),
  TXw ~ theta * Yw,
  TX ~ TXc + TXw,
  Vc ~ Vc[-1] + (YDc - Cc),
  Vw ~ Vw[-1] + (YDw - Cw),
  V ~ Vc + Vw,
  Cc ~ alpha1c * YDc + alpha2c * Vc[-1],
  Cw ~ alpha1w * YDw + alpha2w * Vw[-1],
  C ~ Cc + Cw,
  Hhc ~ Vc - Bh,
  Hhw ~ Vw,
  Hh ~ Hhc + Hhw,
  Hhc1 ~ Vc * ((1 - lambda0) - lambda1 * r + lambda2 * ( YDc/Vc )), # EQ 4.6A
  Hhw1 ~ Vw * ((1 - lambda0) - lambda1 * r  + lambda2 * ( YDw/Vw )), # EQ 4.6A
  Hh1 ~ Hhc1 + Hhw1,
  Bh ~ Vc * (lambda0 + lambda1 * r - lambda2 * ( YDc/Vc )),
  Bs ~ Bs[-1] + (G + r[-1] * Bs[-1]) - (TX + r[-1] * Bcb[-1]),
  Hs ~ Hs[-1] + Bcb - Bcb[-1],
  Bcb ~ Bs - Bh
)

pc_ext <- sfcr_set(
  # Exogenous

  r ~ 0.025,
  G ~ 20,

  # Parameters
  zeta ~ 0.5,
  alpha1c ~ 0.6,
  alpha2c ~ 0.4,
  alpha1w ~ 0.6,
  alpha2w ~ 0.5,
  theta ~ 0.2,
  lambda0 ~ 0.635,
  lambda1 ~ 0.05,
  lambda2 ~ 0.01
)

pc <- sfcr_baseline(
  equations = pc_eqs,
  external = pc_ext,
  periods = 70,
  hidden = c("Hh" = "Hs")
)

sfcr_validate(bs_pc, pc, "bs")

sfcr_validate(tfm_pc, pc, "tfm")

all.equal(pc$Hh, pc$Hh1)

sfcr_dag_cycles_plot(pc_eqs)
sfcr_dag_blocks_plot(pc_eqs)


sfcr_sankey(tfm_pc, pc)


shock1 <- sfcr_shock(
  variables = sfcr_set(
    r ~ 0.035
  ),
  start = 5,
  end = 60
)

pc2 <- sfcr_scenario(pc, scenario = shock1, periods = 60)

pc2 <- pc2 %>%
  mutate(BV = Bh / V,
         McV = Hhc / V,
         MwV = Hhw / V)

pc2_long <- pc2 %>%
  pivot_longer(cols = -period)

pc2_long %>%
  filter(name %in% c("BV", "McV", "MwV")) %>%
  ggplot(aes(x = period, y = value)) +
  geom_line() +
  facet_wrap(~ name, scales = 'free_y') +
  labs(title = "Wealth alocation")

pc2_long %>%
  filter(name %in% c("YDc", "Cc", "YDw", "Cw" )) %>%
  ggplot(aes(x = period, y = value)) +
  geom_line(aes(linetype = name, color = name)) +
  labs(title = "Evolution of disposable income and household consumption",
       subtitle = "Following an increase of 100 points in the rate of interest")



