PROGRAM Solution
IMPLICIT NONE
INTEGER, PARAMETER :: realkind = 8   ! 4: single precision, 8: double precision
REAL(KIND=realkind), PARAMETER :: pi = 3.14159265358979323846264338327950288419716939937510_realkind

REAL(KIND=realkind) :: Q_v, Q_i, Q_2i, D_v, D_i, v_I, v_2i, v_D, v_N, K_iv, r_iv, k_gb, a_0, Omega, b_v
REAL(KIND=realkind) :: f_1, f_2, d_N, P_sl, z_vN, z_iN, z_vL, z_iL

INTEGER :: n_time_step, i_time_step
INTEGER :: n_write_step, i_write_step, i_write
REAL(KIND=realkind) :: dt, dt_half, dt_sixth

REAL(KIND=realkind) :: C_v_factor1, C_v_factor2
REAL(KIND=realkind) :: C_i_factor1, C_i_factor2, C_i_factor3
REAL(KIND=realkind) :: C_2i_factor1, C_2i_factor2, C_2i_factor3
REAL(KIND=realkind) :: rho_N_factor1, rho_N_factor2
REAL(KIND=realkind) :: C_J_factor1, C_J_factor2, C_J_factor3, C_J_factor4
REAL(KIND=realkind) :: d_J_factor1, d_J_factor2, d_J_factor3, d_J_factor4

REAL(KIND=realkind) :: k1_C_v,   k2_C_v,   k3_C_v,   k4_C_v
REAL(KIND=realkind) :: k1_C_i,   k2_C_i,   k3_C_i,   k4_C_i
REAL(KIND=realkind) :: k1_C_2i,  k2_C_2i,  k3_C_2i,  k4_C_2i
REAL(KIND=realkind) :: k1_rho_N, k2_rho_N, k3_rho_N, k4_rho_N
REAL(KIND=realkind) :: k1_C_J,   k2_C_J,   k3_C_J,   k4_C_J
REAL(KIND=realkind) :: k1_d_J,   k2_d_J,   k3_d_J,   k4_d_J

REAL(KIND=realkind) :: t
REAL(KIND=realkind) :: C_v, C_i, C_2i, rho_N, C_J, d_J
REAL(KIND=realkind) :: C_v_ini, C_i_ini, C_2i_ini, rho_N_ini, C_J_ini, d_J_ini
REAL(KIND=realkind) :: C_v_tmp, C_i_tmp, C_2i_tmp, rho_N_tmp, C_J_tmp, d_J_tmp


!!!===   CONSTANTS   ===============================================================================================================
Q_v   = 1.0_realkind
Q_i   = 1.0_realkind
Q_2i  = 1.0_realkind

D_v   = 1.0_realkind
D_i   = 1.0_realkind

v_I   = 1.0_realkind
v_2i  = 1.0_realkind
v_D   = 1.0_realkind
v_N   = 1.0_realkind

K_iv  = 1.0_realkind
r_iv  = 1.0_realkind

k_gb  = 1.0_realkind
a_0   = 1.0_realkind
Omega = 1.0_realkind

b_v   = 1.0_realkind

f_1   = 1.0_realkind
f_2   = 1.0_realkind
d_N   = 1.0_realkind

P_sl  = 1.0_realkind

z_vN  = 1.0_realkind
z_iN  = 1.0_realkind
z_vL  = 1.0_realkind
z_iL  = 1.0_realkind

C_v_ini   = 1.0_realkind
C_i_ini   = 1.0_realkind
C_2i_ini  = 1.0_realkind
rho_N_ini = 1.0_realkind
C_J_ini   = 1.0_realkind
d_J_ini   = 1.0_realkind

t  = 0.0_realkind      ! initial time
dt = 1.0E-3_realkind   ! time increment
n_time_step  = 3000    ! maximum number of time steps
n_write_step = 10      ! write results after every 'n_write_step' time steps

!!!===   INITIALIZATION   ==========================================================================================================
dt_half  = dt / 2
dt_sixth = dt / 6

C_v_factor1 = k_gb ** 2 * D_v
C_v_factor2 = 4 * pi * a_0 ** 2 * D_v / Omega

C_i_factor1 = k_gb ** 2 * D_i
C_i_factor2 = 16 * pi * r_iv * D_i
C_i_factor3 = 4 * pi * a_0 ** 2 * D_v / Omega

C_2i_factor1 =  8 * pi * r_iv * D_i / Omega ** 2
C_2i_factor2 =  4 * pi * r_iv * D_i / Omega
C_2i_factor3 =      pi * v_2i / b_v

rho_N_factor1 = f_1 * ABS(v_D) * pi / d_N
rho_N_factor2 = f_2 * ABS(v_N) / d_N

C_J_factor1 = pi * v_I / (2 * a_0)
C_J_factor2 = 4 * v_I / d_N ** 2
C_J_factor3 = 8 * v_I
C_J_factor4 = f_1 * ABS(v_D) / d_N

d_J_factor1 = 2 * a_0
d_J_factor2 = pi * v_2i / (2 * a_0)
d_J_factor3 = 4 * v_I / (P_sl * d_N ** 2)
d_J_factor4 = 8 * v_I

C_v   = C_v_ini
C_i   = C_i_ini
C_2i  = C_2i_ini
rho_N = rho_N_ini
C_J   = C_J_ini
d_J   = d_J_ini

i_write_step = 0
i_write      = 1
!                123456 23456789 23456789012345678901234 23 23456789012345678901234 23456789012345678901234 23456789012345678901234 23456789012345678901234 23456789012345678901234 23456789012345678901234
WRITE(*, "(A)") "#...i    step           time                          C_v                     C_i                    C_2i&
&                   rho_N                     C_I                     d_J"
WRITE(*, "(I6,I9,ES24.15E3,3x,6ES24.15E3)") i_write, 0, t, C_v, C_i, C_2i, rho_N, C_J, d_J

!!!===   CALCULATIONS   ============================================================================================================
DO i_time_step = 1, n_time_step
   k1_C_v   = Eqn_C_v  (C_v, C_i, C_2i, rho_N, C_J, d_J)
   k1_C_i   = Eqn_C_i  (C_v, C_i, C_2i, rho_N, C_J, d_J)
   k1_C_2i  = Eqn_C_2i (C_v, C_i, C_2i, rho_N, C_J, d_J)
   k1_rho_N = Eqn_rho_N(C_v, C_i, C_2i, rho_N, C_J, d_J)
   k1_C_J   = Eqn_C_J  (C_v, C_i, C_2i, rho_N, C_J, d_J)
   k1_d_J   = Eqn_d_J  (C_v, C_i, C_2i, rho_N, C_J, d_J)

   C_v_tmp   = C_v   + dt_half * k1_C_v
   C_i_tmp   = C_i   + dt_half * k1_C_i
   C_2i_tmp  = C_2i  + dt_half * k1_C_2i
   rho_N_tmp = rho_N + dt_half * k1_rho_N
   C_J_tmp   = C_J   + dt_half * k1_C_J
   d_J_tmp   = d_J   + dt_half * k1_d_J

   k2_C_v   = Eqn_C_v  (C_v_tmp, C_i_tmp, C_2i_tmp, rho_N_tmp, C_J_tmp, d_J_tmp)
   k2_C_i   = Eqn_C_i  (C_v_tmp, C_i_tmp, C_2i_tmp, rho_N_tmp, C_J_tmp, d_J_tmp)
   k2_C_2i  = Eqn_C_2i (C_v_tmp, C_i_tmp, C_2i_tmp, rho_N_tmp, C_J_tmp, d_J_tmp)
   k2_rho_N = Eqn_rho_N(C_v_tmp, C_i_tmp, C_2i_tmp, rho_N_tmp, C_J_tmp, d_J_tmp)
   k2_C_J   = Eqn_C_J  (C_v_tmp, C_i_tmp, C_2i_tmp, rho_N_tmp, C_J_tmp, d_J_tmp)
   k2_d_J   = Eqn_d_J  (C_v_tmp, C_i_tmp, C_2i_tmp, rho_N_tmp, C_J_tmp, d_J_tmp)

   C_v_tmp   = C_v   + dt_half * k2_C_v
   C_i_tmp   = C_i   + dt_half * k2_C_i
   C_2i_tmp  = C_2i  + dt_half * k2_C_2i
   rho_N_tmp = rho_N + dt_half * k2_rho_N
   C_J_tmp   = C_J   + dt_half * k2_C_J
   d_J_tmp   = d_J   + dt_half * k2_d_J

   k3_C_v   = Eqn_C_v  (C_v_tmp, C_i_tmp, C_2i_tmp, rho_N_tmp, C_J_tmp, d_J_tmp)
   k3_C_i   = Eqn_C_i  (C_v_tmp, C_i_tmp, C_2i_tmp, rho_N_tmp, C_J_tmp, d_J_tmp)
   k3_C_2i  = Eqn_C_2i (C_v_tmp, C_i_tmp, C_2i_tmp, rho_N_tmp, C_J_tmp, d_J_tmp)
   k3_rho_N = Eqn_rho_N(C_v_tmp, C_i_tmp, C_2i_tmp, rho_N_tmp, C_J_tmp, d_J_tmp)
   k3_C_J   = Eqn_C_J  (C_v_tmp, C_i_tmp, C_2i_tmp, rho_N_tmp, C_J_tmp, d_J_tmp)
   k3_d_J   = Eqn_d_J  (C_v_tmp, C_i_tmp, C_2i_tmp, rho_N_tmp, C_J_tmp, d_J_tmp)

   C_v_tmp   = C_v   + dt * k3_C_v
   C_i_tmp   = C_i   + dt * k3_C_i
   C_2i_tmp  = C_2i  + dt * k3_C_2i
   rho_N_tmp = rho_N + dt * k3_rho_N
   C_J_tmp   = C_J   + dt * k3_C_J
   d_J_tmp   = d_J   + dt * k3_d_J

   k4_C_v   = Eqn_C_v  (C_v_tmp, C_i_tmp, C_2i_tmp, rho_N_tmp, C_J_tmp, d_J_tmp)
   k4_C_i   = Eqn_C_i  (C_v_tmp, C_i_tmp, C_2i_tmp, rho_N_tmp, C_J_tmp, d_J_tmp)
   k4_C_2i  = Eqn_C_2i (C_v_tmp, C_i_tmp, C_2i_tmp, rho_N_tmp, C_J_tmp, d_J_tmp)
   k4_rho_N = Eqn_rho_N(C_v_tmp, C_i_tmp, C_2i_tmp, rho_N_tmp, C_J_tmp, d_J_tmp)
   k4_C_J   = Eqn_C_J  (C_v_tmp, C_i_tmp, C_2i_tmp, rho_N_tmp, C_J_tmp, d_J_tmp)
   k4_d_J   = Eqn_d_J  (C_v_tmp, C_i_tmp, C_2i_tmp, rho_N_tmp, C_J_tmp, d_J_tmp)

   t = t + dt
   C_v   = C_v   + (k1_C_v   + 2 * (k2_C_v   + k3_C_v  ) + k4_C_v  ) * dt_sixth
   C_i   = C_i   + (k1_C_i   + 2 * (k2_C_i   + k3_C_i  ) + k4_C_i  ) * dt_sixth
   C_2i  = C_2i  + (k1_C_2i  + 2 * (k2_C_2i  + k3_C_2i ) + k4_C_2i ) * dt_sixth
   rho_N = rho_N + (k1_rho_N + 2 * (k2_rho_N + k3_rho_N) + k4_rho_N) * dt_sixth
   C_J   = C_J   + (k1_C_J   + 2 * (k2_C_J   + k3_C_J  ) + k4_C_J  ) * dt_sixth
   d_J   = d_J   + (k1_d_J   + 2 * (k2_d_J   + k3_d_J  ) + k4_d_J  ) * dt_sixth

   IF (n_write_step > 0) THEN
      i_write_step = i_write_step + 1
      IF (i_write_step == n_write_step) THEN
         i_write = i_write + 1
         WRITE(*, "(I6,I9,ES24.15E3,3x,6ES24.15E3)") i_write, i_time_step, t, C_v, C_i, C_2i, rho_N, C_J, d_J
         i_write_step = 0
      END IF
   END IF
END DO

IF ((n_write_step <= 0) .OR. (i_write_step > 0)) THEN
   i_write = i_write + 1
   WRITE(*, "(I6,I9,ES24.15E3,3x,6ES24.15E3)") i_write, n_time_step, t, C_v, C_i, C_2i, rho_N, C_J, d_J
END IF


CONTAINS

!!!#################################################################################################################################
REAL(KIND=realkind) FUNCTION K_vN(x)
REAL(KIND=realkind), INTENT(IN) :: x

K_vN = z_vN * x

END FUNCTION K_vN

!!!#################################################################################################################################
REAL(KIND=realkind) FUNCTION K_iN(x)
REAL(KIND=realkind), INTENT(IN) :: x

K_iN = z_iN * x

END FUNCTION K_iN

!!!#################################################################################################################################
REAL(KIND=realkind) FUNCTION K_vL(x)
REAL(KIND=realkind), INTENT(IN) :: x

K_vL = z_vL * x

END FUNCTION K_vL

!!!#################################################################################################################################
REAL(KIND=realkind) FUNCTION K_iL(x)
REAL(KIND=realkind), INTENT(IN) :: x

K_iL = z_iL * x

END FUNCTION K_iL

!!!#################################################################################################################################
REAL(KIND=realkind) FUNCTION Eqn_C_v(C_v, C_i, C_2i, rho_N, C_J, d_J)
REAL(KIND=realkind), INTENT(IN) :: C_v, C_i, C_2i, rho_N, C_J, d_J

Eqn_C_v = Q_v - (K_iv * C_i + K_vN(rho_N) * D_v + K_vL(pi * C_J * d_J) * D_v + C_v_factor1 + C_v_factor2 * C_2i) * C_v

END FUNCTION Eqn_C_v

!!!#################################################################################################################################
REAL(KIND=realkind) FUNCTION Eqn_C_i(C_v, C_i, C_2i, rho_N, C_J, d_J)
REAL(KIND=realkind), INTENT(IN) :: C_v, C_i, C_2i, rho_N, C_J, d_J

Eqn_C_i = Q_i - (K_iv * C_v + K_iN(rho_N) * D_i + K_iL(pi * C_J * d_J) * D_i + C_i_factor1 + C_i_factor2 * C_i) * C_i &
        + C_i_factor3 * C_v * C_2i
!

END FUNCTION Eqn_C_i

!!!#################################################################################################################################
REAL(KIND=realkind) FUNCTION Eqn_C_2i(C_v, C_i, C_2i, rho_N, C_J, d_J)
REAL(KIND=realkind), INTENT(IN) :: C_v, C_i, C_2i, rho_N, C_J, d_J

Eqn_C_2i = Q_2i + C_2i_factor1 * C_i ** 2 - (C_2i_factor2 * C_v + C_2i_factor3) * C_2i

END FUNCTION Eqn_C_2i

!!!#################################################################################################################################
REAL(KIND=realkind) FUNCTION Eqn_rho_N(C_v, C_i, C_2i, rho_N, C_J, d_J)
REAL(KIND=realkind), INTENT(IN) :: C_v, C_i, C_2i, rho_N, C_J, d_J

Eqn_rho_N = rho_N_factor1 * d_J * C_J - rho_N_factor2 * rho_N

END FUNCTION Eqn_rho_N

!!!#################################################################################################################################
REAL(KIND=realkind) FUNCTION Eqn_C_J(C_v, C_i, C_2i, rho_N, C_J, d_J)
REAL(KIND=realkind), INTENT(IN) :: C_v, C_i, C_2i, rho_N, C_J, d_J

Eqn_C_J = C_J_factor1 * C_2i - ((C_J_factor2 + C_J_factor3 * d_J * C_J) * d_J + C_J_factor4) * C_J

END FUNCTION Eqn_C_J

!!!#################################################################################################################################
REAL(KIND=realkind) FUNCTION Eqn_d_J(C_v, C_i, C_2i, rho_N, C_J, d_J)
REAL(KIND=realkind), INTENT(IN) :: C_v, C_i, C_2i, rho_N, C_J, d_J

Eqn_d_J = v_I - (d_J - d_J_factor1) * d_J_factor2 * C_2i / C_J - ((d_N - d_J) * d_J_factor3 + d_J_factor4 * d_J ** 2 * C_J) * d_J

END FUNCTION Eqn_d_J

!!!#################################################################################################################################

END PROGRAM Solution


