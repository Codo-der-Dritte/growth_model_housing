bmw_eqs <- sfcr_set(
  # Basic behavioral equations
  Cs ~ Cd,
  Is ~ Id,
  Ns ~ Nd,
  Ls ~ Ls[-1] + Ld - Ld[-1],

  # Transactions of the firms
  Y ~ Cs + Is,
  WBd ~ Y - rl[-1] * Ld[-1] - AF,
  AF ~ delta * K[-1],
  Ld ~ Ld[-1] + Id - AF,

  # Transactions of households
  YD ~ WBs + rm[-1] * Mh[-1],
  Mh ~ Mh[-1] + YD - Cd,

  # Transactions of the banks
  Ms ~ Ms[-1] + Ls - Ls[-1],
  rm ~ rl,

  # The wage bill
  WBs ~ W * Ns,
  Nd ~ Y / pr,
  W ~ WBd / Nd,

  # Household behavior
  Cd ~ alpha0 + alpha1 * YD + alpha2 * Mh[-1],

  # The investment behavior
  K ~ K[-1] + Id - DA,
  DA ~ delta * K[-1],
  KT ~ kappa * Y[-1],
  Id ~ gamma * (KT - K[-1]) + DA
)

bmw_external <- sfcr_set(
  rl ~ 0.025,
  alpha0 ~ 20,
  alpha1 ~ 0.75,
  alpha2 ~ 0.10,
  delta ~ 0.10,
  gamma ~ 0.15,
  kappa ~ 1,
  pr ~ 1
)
