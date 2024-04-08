
# Install and Load libraries

libs <- c("tidyverse", "leaflet", "lubridate", "readxl", "scales", "xlsx",
          "ggraph", "networkD3", "igraph", "expm", "sfcr")
installed_libs <- libs %in% rownames(installed.packages())
if (any(installed_libs == F)) {
  install.packages(libs[!installed_libs])
}
invisible(lapply(libs, library, character.only = T))

# Start with INSOUT Model of Monetary economics.

bs_insout <- sfcr_matrix(
  columns = c("Households", "Firms", "Government", "Central bank", "Banks", "Sum"),
  codes = c("h", "f", "g", "cb", "b", "s"),
  r1 = c("Inventories", f = "+INV", s = "+INV"),
  r2 = c("HPM", h = "+Hhd", cb = "-Hs", b = "+Hbd"),
  r3 = c("Advances", cb = "+As", b = "-Ad"),
  r4 = c("Checking deposits", h = "+M1h", b = "-M1s"),
  r5 = c("Bills", h = "+Bhh", g = "-Bs", cb = "+Bcb", b = "+Bbd"),
  r6 = c("Loans", f = "-Ld", b = "+Ls"),
  r7 = c("Balance", h = "-V", f = 0, g = "+GD", cb = 0, b = 0, s = "-INV")
)

sfcr_matrix_display(bs_insout, "bs")

tfm_insout <- sfcr_matrix(
  columns = c("Households", "Firms current", "Firms capital", "Govt.", "CB current", "CB capital", "Banks current", "Banks capital"),
  codes = c("h", "fc", "fk", "g", "cbc", "cbk", "bc", "bk"),
  c("Consumption", h = "-C", fc = "+C"),
  c("Govt. Expenditures", fc = "+G", g = "-G"),
  c("Ch. Inv", fc = "+(INV - INV[-1])", fk = "-(INV - INV[-1])"),
  c("Taxes", fc = "-TX", g = "+TX"),
  c("Wages", h = "+WB", fc = "-WB"),
  c("Entrepreneurial profits", h = "+FXf", fc = "-FXf"),
  c("Bank profits", h = "+FXb", bc = "-FXb"),
  c("CB profits", g = "+FXcb", cbc = "-FXcb"),
  c("int. advances", cbc = "+ra[-1] * As[-1]", bc = "-ra[-1] * Ad[-1]"),
  c("int. loans", fc = "-rl[-1] * Ld[-1]", bc = "+rl[-1] * Ld[-1]"),
  c("int. deposits", h = "+rm[-1] * M1h[-1]", bc = "-rm[-1] * M1h[-1]"),
  c("int. bills", h = "+rb[-1] * Bhh[-1]", g = "-rb[-1] * Bs[-1]", cbc = "+rb[-1] * Bcb[-1]", bc = "+rb[-1] * Bbd[-1]"),
  c("Ch. advances", cbk = "-(As - As[-1])", bk = "+(Ad - Ad[-1])"),
  c("Ch. loans", fk = "+(Ld - Ld[-1])", bk = "-(Ls - Ls[-1])"),
  c("Ch. cash", h = "-(Hhh - Hhh[-1])", cbk = "+(Hs - Hs[-1])", bk = "-(Hbd - Hbd[-1])"),
  c("Ch. M1", h = "-(M1h - M1h[-1])", bk = "+(M1s - M1s[-1])"),
  c("Ch. bills", h = "-(Bhh - Bhh[-1])", g = "+(Bs - Bs[-1])", cbk = "-(Bcb - Bcb[-1])", bk = "-(Bbd - Bbd[-1])"),
)
sfcr_matrix_display(tfm_insout, "tfm")


insout_eqs <- sfcr_set(

  # Firm's behavioral equations

  y ~ sE + (invE - inv[-1]),
  N ~ y / pr,
  WB ~ N * W,
  UC ~ WB / y,
  sE ~ beta * s[-1] + (1 - beta) * sE[-1],
  invT ~ sigmaT * sE,
  sigmaT ~ sigma0 - sigma1 * rl,
  #rrl ~ ((1 + rl) / (1 + pi)) - 1,
  invE ~ inv[-1] + gamma * (invT - inv[-1]),
  p ~ (1 + tau) * (1 + phi) * NHUC,
  NHUC ~ (1 - sigmaT) * UC + sigmaT * (1 + rl[-1]) * UC[-1], ## rl[-1] instead of rl
  # FXfE ~ (phi / (1 + phi)) * (1 / (1 + tau)) * p * sE,
  # OK

  # Firm's realized outcomes

  s ~ c + g,
  S ~ s * p,
  inv ~ inv[-1] + y - s,
  sigmas ~ inv[-1] / s,
  INV ~ inv * UC,
  Ld ~ INV,
  FXf ~ S - TX - WB + (INV - INV[-1]) - rl[-1] * INV[-1], # rl[-1] instead of rl
  pi ~ (p / p[-1]) - 1,
  # OK

  # Households realized outcomes

  YDr ~ FX + WB + rm[-1] * M1h[-1] + rb[-1] * Bhh[-1], ## Here's a mistake in Zezza's code. It should be M2h[-1] as in G&L and NOT M2d.
  #CG ~ (pbl - pbl[-1]) * BLh[-1],
  #YDhs ~ YDr + CG,
  FX ~ FXf + FXb,
  V ~ V[-1] + YDhs - C,
  Vnc ~ V - Hhh, ## Zezza writes it as Hhd instead of Hhh. Here it is harmless, but theoretically it should be Hhh and not Hhd as what matters for wealth net of cash is the realized holdings of cash.
  ydr ~ YDr/p - pi * (V[-1]/p),
  #ydhs ~ (YDr - pi * V[-1] + CG) / p,
  #ydhs ~ c + v - v[-1], # Equation 10.27A
  v ~ V/p,
  # OK

  # Households behavioral

  c ~ alpha0 + alpha1 * ydrE + alpha2 * v[-1],
  ydrE ~ epsilon * ydr[-1] + (1 - epsilon) * ydrE[-1],
  C ~ p * c,
  YDrE ~ p * ydrE + pi * (V[-1]/p),
  VE ~ V[-1] + (YDrE - C),
  Hhd ~ lambdac * C,
  VncE ~ VE - Hhd,
  ERrbl ~ rbl, ## ERrbl is not on the list of equations. I kept it simple.
  # OK

  # Households' portfolio equations

  # There's no M1d equation in Zezza's code as they are not necessary
  M2d ~ VncE * (lambda20 + lambda22 * rm + lambda23 * rb + lambda24 * ERrbl + lambda25 * (YDrE / VncE)),
  Bhd ~ VncE * (lambda30 + lambda32 * rm + lambda33 * rb + lambda34 * ERrbl + lambda35 * (YDrE / VncE)),
  BLd ~ (VncE / pbl) * (lambda40 + lambda42 * rm + lambda43 * rb + lambda44 * ERrbl + lambda45 * (YDrE / VncE)),
  # OK

  # However, it is a good exercise to write them down. If all values are correct, the equations
  # below must be equal.
  M1d ~ VncE * (lambda10 + lambda12 * rm + lambda13 * rb + lambda14 * ERrbl + lambda15 * (YDrE / VncE)),
  M1d2 ~ VncE - M2d - Bhd - pbl * BLd,


  # Realized portfolio asset holdings

  # Bhs ~ Bhd, # Not explicit in GL
  Hhh ~ Hhd,
  Bhh ~ Bhd,
  #BLh ~ BLd,
  M1hN ~ Vnc - M2d - Bhd - pbl * BLd,
  z1 ~ if (M1hN > 0) {1} else {0},
  z2 ~ 1 - z1,
  M1h ~ M1hN * z1,
  M2hN ~ M2d,
  M2h ~ M2d * z1 + (Vnc - Bhh - pbl * BLd) * z2,
  # OK

  # Government's equations

  TX ~ S * (tau / (1 + tau)),
  G ~ p * g,
  PSBR ~ G + rb[-1] * Bs[-1] + BLs[-1] - (TX + FXcb),
  Bs ~ Bs[-1] + PSBR - (BLs - BLs[-1]) * pbl,
  BLs ~ BLd,
  pbl ~ 1 / rbl,

  GD ~ GD[-1] + PSBR, # Not in the equations list, but I added to check BS consistency
  # OK

  # Central bank's equations

  Hs ~ Bcb + As,
  Hbs ~ Hs - Hhs,
  Bcb ~ Bs - Bhh - Bbd,
  As ~ Ad,
  ra ~ rb,
  FXcb ~ rb[-1] * Bcb[-1] + ra[-1] * As[-1],
  # OK

  # Bank's realized (supply) equations

  Hhs ~ Hhd,
  M1s ~ M1h, ## M1h instead of M1d as in Zezza. Bank's supply the realized portfolio holdings and NOT its notional demand.
  M2s ~ M2h, ## M2h instead of M2d as in Zezza. Bank's supply the realized portfolio holdings and NOT its notional demand.

  # Otherwise, the model would not close because M1d depends on EXPECTED wealth but the REALIZED demand for deposits depends
  # on REALIZED wealth. It is a residual variable. The bank's supply the deposits that are ACTUALLY needed.
  # You can check that M1d, however measured, only equals M1h in the stationary state.

  Ls ~ Ld,
  Hbd ~ ro1 * M1s + ro2 * M2s,
  # OK

  # Bank's balance sheet constraints

  BbdN ~ M1s + M2s - Ls - Hbd,
  BLRN ~ BbdN / (M1s + M2s),
  Ad ~ (bot * (M1s + M2s) - BbdN) * z3, ## Z3 instead of Z4
  z3 ~ if (BLRN < bot) {1} else {0},
  Bbd ~ Ad + M1s + M2s - Ls - Hbd,
  BLR ~ Bbd / (M1s + M2s),
  # OK

  # Determination of interest rates by banks

  rm ~ rm[-1] + zetam * (z4 - z5) + zetab * (rb - rb[-1]),
  z4 ~ if (BLRN[-1] < bot) {1} else {0},
  z5 ~ if (BLRN[-1] > top) {1} else {0},
  FXb ~ rl[-1] * Ls[-1] + rb[-1] * Bbd[-1] - rm[-1] * M2s[-1] - ra[-1] * Ad[-1],
  rl ~ rl[-1] + zetal * (z6 - z7) + (rb - rb[-1]),
  z6 ~ if (BPM < botpm) {1} else {0},
  z7 ~ if (BPM > toppm) {1} else {0},
  # Since the sfcr package does not accept more than 1 lag directly,
  # we need to create two auxiliary values that are lags of endogenous
  # variables and then take the lag of these variables
  lM1s ~ M1s[-1],
  lM2s ~ M2s[-1],
  BPM ~ (FXb + FXb[-1]) / (lM1s + lM1s[-1] + lM2s + lM2s[-1]),
  # OK

  # Inflationary forces

  #omegaT ~ Omega0 + Omega1 * pr + Omega2 * (N / Nfe),
  omegaT ~ exp(Omega0 + Omega1 * log(pr) + Omega2 * log(N / Nfe)), ## Zezza's equation

  # The model explodes if G&L's definition of omegaT is used.

  W ~ W[-1] * (1 + Omega3 * (omegaT[-1] - (W[-1] / p[-1]))),
  # Nfe ~ s / pr, ## Zezza's definition in model DISINF. Doesn't work here.
  # Nfe ~ s[-1] / pr, ## One possible solution
  Y ~ p * s + UC * (inv - inv[-1])

)

insout_ext <- sfcr_set(
  # EXOGENOUS

  rbl ~ 0.027,
  rb ~ 0.023,
  pr ~ 1,
  g ~ 25,
  Nfe ~ 133.28, ## Zezza supplies this exogenous full employment value
  #Nfe ~ 142.75,

  # In model DISINF Nfe was defined as `s / pr`, but it doesn't work here.
  # I discuss this issue in the scenarios.

  # PARAMETERS

  alpha0 ~ 0,
  alpha1 ~ 0.95,
  alpha2 ~ 0.05,
  beta ~ 0.5,
  bot = bot ~ 0.02,
  botpm = botpm ~ 0.003,
  epsilon ~ 0.5,
  gamma ~ 0.5,
  lambdac ~ 0.1,
  phi ~ 0.1,
  ro1 ~ 0.1,
  ro2 ~ 0.1,
  sigma0 ~ 0.3612,
  sigma1 ~ 3,
  tau ~ 0.25,
  zetab ~  0.9,
  zetal ~  0.0002,
  zetam ~  0.0002,
  Omega0 ~ -0.32549,
  Omega1 ~ 1,
  Omega2 ~ 1.5,
  Omega3 ~ 0.1,
  #top ~ 0.04, # Zezza's top
  top = top ~  0.06,
  toppm = toppm ~ 0.005,



  # PORTFOLIO PARAMETERS

  # (note that the sfcr package does not evaluate variables defined in the global environment when setting the external variables.
  # Therefore, the parameters' values must be added by hand here.)

  lambda10 ~ -0.17071,
  lambda11 ~ 0,
  lambda12 ~ 0,
  lambda13 ~ 0,
  lambda14 ~ 0,
  lambda15 ~ 0.18,
  lambda20 ~ 0.52245,
  lambda21 ~ 0,
  lambda22 ~ 30,
  lambda23 ~ -15,
  lambda24 ~ -15,
  lambda25 ~ -0.06,
  lambda30 ~ 0.47311,
  lambda31 ~ 0,
  lambda32 ~ -15,
  lambda33 ~ 30,
  lambda34 ~ -15,
  lambda35 ~ -0.06,
  lambda40 ~ 0.17515,
  lambda41 ~ 0,
  lambda42 ~ -15,
  lambda43 ~ -15,
  lambda44 ~ 30,
  lambda45 ~ -0.06

)

insout <- sfcr_baseline(
  equations = insout_eqs,
  external = insout_ext,
  periods = 210,
  initial = sfcr_set(p ~ 1, W ~ 1, UC ~ 1, BPM ~ 0.0035),
  hidden = c("Hbd" = "Hbs"),
  tol = 1e-15
)

insout %>%
  pivot_longer(cols = -period) %>%
  filter(name %in% c("Y", "y", "s", "inv", "pi", "Bs", "M1s", "M2s", "V", "INV", "FXf", "FXb")) %>%
  ggplot(aes(x = period, y = value)) +
  geom_line() +
  facet_wrap(~name, scales = "free_y") +
  labs(title = "INSOUT", subtitle = "Stable steady state")

sfcr_validate(bs_insout, insout, which = "bs")
sfcr_validate(tfm_insout, insout, which = "tfm")









