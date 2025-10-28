#lang racket

(require "../prover_core.rkt")
(require "../prover_vector.rkt")

 (prove-roe-vector2-1d-flux-conservation '#hash((cons-exprs . (mom_y mom_z)) (flux-exprs . ((* mom_y u) (* mom_z u))) (max-speed-exprs . ((abs u) (abs u))) (name . "isothermal-euler-mom-yz") (parameters . ((define u 0.0)))) #:cfl 0.95 #:init-funcs '(0.0 0.0) #:nx 200 #:t-final 0.1 #:x0 0.0 #:x1 1.0)
  (symbolic-jacobian '((* mom_y u) (* mom_z u)) '(mom_y mom_z))
   (symbolic-diff '(* mom_y u) 'mom_y)
    (symbolic-diff 'mom_y 'mom_y)
    1.0
    (symbolic-diff 'u 'mom_y)
    0.0
   '(+ (* 1.0 u) (* mom_y 0.0))
   (symbolic-simp '(+ (* 1.0 u) (* mom_y 0.0)))
    (symbolic-simp-rule '(+ (* 1.0 u) (* mom_y 0.0)))
     (symbolic-simp-rule '(* 1.0 u))
     'u
     (symbolic-simp-rule '(* mom_y 0.0))
     '(* 0.0 mom_y)
    '(+ u (* 0.0 mom_y))
   (symbolic-simp '(+ u (* 0.0 mom_y)))
    (symbolic-simp-rule '(+ u (* 0.0 mom_y)))
     (symbolic-simp-rule 'u)
     'u
     (symbolic-simp-rule '(* 0.0 mom_y))
     0.0
    '(+ u 0.0)
   (symbolic-simp '(+ u 0.0))
    (symbolic-simp-rule '(+ u 0.0))
    '(+ 0.0 u)
   (symbolic-simp '(+ 0.0 u))
    (symbolic-simp-rule '(+ 0.0 u))
    'u
   (symbolic-simp 'u)
    (symbolic-simp-rule 'u)
    'u
   'u
   (symbolic-diff '(* mom_y u) 'mom_z)
    (symbolic-diff 'mom_y 'mom_z)
    0.0
    (symbolic-diff 'u 'mom_z)
    0.0
   '(+ (* 0.0 u) (* mom_y 0.0))
   (symbolic-simp '(+ (* 0.0 u) (* mom_y 0.0)))
    (symbolic-simp-rule '(+ (* 0.0 u) (* mom_y 0.0)))
     (symbolic-simp-rule '(* 0.0 u))
     0.0
     (symbolic-simp-rule '(* mom_y 0.0))
     '(* 0.0 mom_y)
    '(+ 0.0 (* 0.0 mom_y))
   (symbolic-simp '(+ 0.0 (* 0.0 mom_y)))
    (symbolic-simp-rule '(+ 0.0 (* 0.0 mom_y)))
    '(* 0.0 mom_y)
   (symbolic-simp '(* 0.0 mom_y))
    (symbolic-simp-rule '(* 0.0 mom_y))
    0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff '(* mom_z u) 'mom_y)
    (symbolic-diff 'mom_z 'mom_y)
    0.0
    (symbolic-diff 'u 'mom_y)
    0.0
   '(+ (* 0.0 u) (* mom_z 0.0))
   (symbolic-simp '(+ (* 0.0 u) (* mom_z 0.0)))
    (symbolic-simp-rule '(+ (* 0.0 u) (* mom_z 0.0)))
     (symbolic-simp-rule '(* 0.0 u))
     0.0
     (symbolic-simp-rule '(* mom_z 0.0))
     '(* 0.0 mom_z)
    '(+ 0.0 (* 0.0 mom_z))
   (symbolic-simp '(+ 0.0 (* 0.0 mom_z)))
    (symbolic-simp-rule '(+ 0.0 (* 0.0 mom_z)))
    '(* 0.0 mom_z)
   (symbolic-simp '(* 0.0 mom_z))
    (symbolic-simp-rule '(* 0.0 mom_z))
    0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff '(* mom_z u) 'mom_z)
    (symbolic-diff 'mom_z 'mom_z)
    1.0
    (symbolic-diff 'u 'mom_z)
    0.0
   '(+ (* 1.0 u) (* mom_z 0.0))
   (symbolic-simp '(+ (* 1.0 u) (* mom_z 0.0)))
    (symbolic-simp-rule '(+ (* 1.0 u) (* mom_z 0.0)))
     (symbolic-simp-rule '(* 1.0 u))
     'u
     (symbolic-simp-rule '(* mom_z 0.0))
     '(* 0.0 mom_z)
    '(+ u (* 0.0 mom_z))
   (symbolic-simp '(+ u (* 0.0 mom_z)))
    (symbolic-simp-rule '(+ u (* 0.0 mom_z)))
     (symbolic-simp-rule 'u)
     'u
     (symbolic-simp-rule '(* 0.0 mom_z))
     0.0
    '(+ u 0.0)
   (symbolic-simp '(+ u 0.0))
    (symbolic-simp-rule '(+ u 0.0))
    '(+ 0.0 u)
   (symbolic-simp '(+ 0.0 u))
    (symbolic-simp-rule '(+ 0.0 u))
    'u
   (symbolic-simp 'u)
    (symbolic-simp-rule 'u)
    'u
   'u
  '((u 0.0) (0.0 u))
  (symbolic-roe-matrix '((u 0.0) (0.0 u)) '(mom_y mom_z))
   (flux-deriv-replace 'u 'mom_y 'mom_yL)
   'u
   (flux-deriv-replace 'u 'mom_z 'mom_zL)
   'u
   (flux-deriv-replace 'u 'mom_y 'mom_yR)
   'u
   (flux-deriv-replace 'u 'mom_z 'mom_zR)
   'u
   (symbolic-simp '(+ (* 0.5 u) (* 0.5 u)))
    (symbolic-simp-rule '(+ (* 0.5 u) (* 0.5 u)))
    '(* (+ 0.5 0.5) u)
   (symbolic-simp '(* (+ 0.5 0.5) u))
    (symbolic-simp-rule '(* (+ 0.5 0.5) u))
     (symbolic-simp-rule '(+ 0.5 0.5))
     1.0
     (symbolic-simp-rule 'u)
     'u
    '(* 1.0 u)
   (symbolic-simp '(* 1.0 u))
    (symbolic-simp-rule '(* 1.0 u))
    'u
   (symbolic-simp 'u)
    (symbolic-simp-rule 'u)
    'u
   'u
   (flux-deriv-replace 0.0 'mom_y 'mom_yL)
   0.0
   (flux-deriv-replace 0.0 'mom_z 'mom_zL)
   0.0
   (flux-deriv-replace 0.0 'mom_y 'mom_yR)
   0.0
   (flux-deriv-replace 0.0 'mom_z 'mom_zR)
   0.0
   (symbolic-simp '(+ (* 0.5 0.0) (* 0.5 0.0)))
    (symbolic-simp-rule '(+ (* 0.5 0.0) (* 0.5 0.0)))
    '(* (+ 0.5 0.5) 0.0)
   (symbolic-simp '(* (+ 0.5 0.5) 0.0))
    (symbolic-simp-rule '(* (+ 0.5 0.5) 0.0))
    '(* 0.0 (+ 0.5 0.5))
   (symbolic-simp '(* 0.0 (+ 0.5 0.5)))
    (symbolic-simp-rule '(* 0.0 (+ 0.5 0.5)))
    0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (flux-deriv-replace 0.0 'mom_y 'mom_yL)
   0.0
   (flux-deriv-replace 0.0 'mom_z 'mom_zL)
   0.0
   (flux-deriv-replace 0.0 'mom_y 'mom_yR)
   0.0
   (flux-deriv-replace 0.0 'mom_z 'mom_zR)
   0.0
   (symbolic-simp '(+ (* 0.5 0.0) (* 0.5 0.0)))
    (symbolic-simp-rule '(+ (* 0.5 0.0) (* 0.5 0.0)))
    '(* (+ 0.5 0.5) 0.0)
   (symbolic-simp '(* (+ 0.5 0.5) 0.0))
    (symbolic-simp-rule '(* (+ 0.5 0.5) 0.0))
    '(* 0.0 (+ 0.5 0.5))
   (symbolic-simp '(* 0.0 (+ 0.5 0.5)))
    (symbolic-simp-rule '(* 0.0 (+ 0.5 0.5)))
    0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (flux-deriv-replace 'u 'mom_y 'mom_yL)
   'u
   (flux-deriv-replace 'u 'mom_z 'mom_zL)
   'u
   (flux-deriv-replace 'u 'mom_y 'mom_yR)
   'u
   (flux-deriv-replace 'u 'mom_z 'mom_zR)
   'u
   (symbolic-simp '(+ (* 0.5 u) (* 0.5 u)))
    (symbolic-simp-rule '(+ (* 0.5 u) (* 0.5 u)))
    '(* (+ 0.5 0.5) u)
   (symbolic-simp '(* (+ 0.5 0.5) u))
    (symbolic-simp-rule '(* (+ 0.5 0.5) u))
     (symbolic-simp-rule '(+ 0.5 0.5))
     1.0
     (symbolic-simp-rule 'u)
     'u
    '(* 1.0 u)
   (symbolic-simp '(* 1.0 u))
    (symbolic-simp-rule '(* 1.0 u))
    'u
   (symbolic-simp 'u)
    (symbolic-simp-rule 'u)
    'u
   'u
  '((u 0.0) (0.0 u))
  (symbolic-simp '(- mom_yL mom_yR))
   (symbolic-simp-rule '(- mom_yL mom_yR))
    (symbolic-simp-rule 'mom_yL)
    'mom_yL
    (symbolic-simp-rule 'mom_yR)
    'mom_yR
   '(- mom_yL mom_yR)
  '(- mom_yL mom_yR)
  (symbolic-simp '(- mom_zL mom_zR))
   (symbolic-simp-rule '(- mom_zL mom_zR))
    (symbolic-simp-rule 'mom_zL)
    'mom_zL
    (symbolic-simp-rule 'mom_zR)
    'mom_zR
   '(- mom_zL mom_zR)
  '(- mom_zL mom_zR)
  (symbolic-simp '(+ (* u (- mom_yL mom_yR)) (* 0.0 (- mom_zL mom_zR))))
   (symbolic-simp-rule '(+ (* u (- mom_yL mom_yR)) (* 0.0 (- mom_zL mom_zR))))
    (symbolic-simp-rule '(* u (- mom_yL mom_yR)))
     (symbolic-simp-rule 'u)
     'u
     (symbolic-simp-rule '(- mom_yL mom_yR))
      (symbolic-simp-rule 'mom_yL)
      'mom_yL
      (symbolic-simp-rule 'mom_yR)
      'mom_yR
     '(- mom_yL mom_yR)
    '(* u (- mom_yL mom_yR))
    (symbolic-simp-rule '(* 0.0 (- mom_zL mom_zR)))
    0.0
   '(+ (* u (- mom_yL mom_yR)) 0.0)
  (symbolic-simp '(+ (* u (- mom_yL mom_yR)) 0.0))
   (symbolic-simp-rule '(+ (* u (- mom_yL mom_yR)) 0.0))
   '(+ 0.0 (* u (- mom_yL mom_yR)))
  (symbolic-simp '(+ 0.0 (* u (- mom_yL mom_yR))))
   (symbolic-simp-rule '(+ 0.0 (* u (- mom_yL mom_yR))))
   '(* u (- mom_yL mom_yR))
  (symbolic-simp '(* u (- mom_yL mom_yR)))
   (symbolic-simp-rule '(* u (- mom_yL mom_yR)))
    (symbolic-simp-rule 'u)
    'u
    (symbolic-simp-rule '(- mom_yL mom_yR))
     (symbolic-simp-rule 'mom_yL)
     'mom_yL
     (symbolic-simp-rule 'mom_yR)
     'mom_yR
    '(- mom_yL mom_yR)
   '(* u (- mom_yL mom_yR))
  '(* u (- mom_yL mom_yR))
  (symbolic-simp '(+ (* 0.0 (- mom_yL mom_yR)) (* u (- mom_zL mom_zR))))
   (symbolic-simp-rule '(+ (* 0.0 (- mom_yL mom_yR)) (* u (- mom_zL mom_zR))))
    (symbolic-simp-rule '(* 0.0 (- mom_yL mom_yR)))
    0.0
    (symbolic-simp-rule '(* u (- mom_zL mom_zR)))
     (symbolic-simp-rule 'u)
     'u
     (symbolic-simp-rule '(- mom_zL mom_zR))
      (symbolic-simp-rule 'mom_zL)
      'mom_zL
      (symbolic-simp-rule 'mom_zR)
      'mom_zR
     '(- mom_zL mom_zR)
    '(* u (- mom_zL mom_zR))
   '(+ 0.0 (* u (- mom_zL mom_zR)))
  (symbolic-simp '(+ 0.0 (* u (- mom_zL mom_zR))))
   (symbolic-simp-rule '(+ 0.0 (* u (- mom_zL mom_zR))))
   '(* u (- mom_zL mom_zR))
  (symbolic-simp '(* u (- mom_zL mom_zR)))
   (symbolic-simp-rule '(* u (- mom_zL mom_zR)))
    (symbolic-simp-rule 'u)
    'u
    (symbolic-simp-rule '(- mom_zL mom_zR))
     (symbolic-simp-rule 'mom_zL)
     'mom_zL
     (symbolic-simp-rule 'mom_zR)
     'mom_zR
    '(- mom_zL mom_zR)
   '(* u (- mom_zL mom_zR))
  '(* u (- mom_zL mom_zR))
  (flux-deriv-replace '(* mom_y u) 'mom_y 'mom_yL)
   (flux-deriv-replace 'mom_y 'mom_y 'mom_yL)
   'mom_yL
   (flux-deriv-replace 'u 'mom_y 'mom_yL)
   'u
  '(* mom_yL u)
  (flux-deriv-replace '(* mom_yL u) 'mom_z 'mom_zL)
   (flux-deriv-replace 'mom_yL 'mom_z 'mom_zL)
   'mom_yL
   (flux-deriv-replace 'u 'mom_z 'mom_zL)
   'u
  '(* mom_yL u)
  (flux-deriv-replace '(* mom_y u) 'mom_y 'mom_yR)
   (flux-deriv-replace 'mom_y 'mom_y 'mom_yR)
   'mom_yR
   (flux-deriv-replace 'u 'mom_y 'mom_yR)
   'u
  '(* mom_yR u)
  (flux-deriv-replace '(* mom_yR u) 'mom_z 'mom_zR)
   (flux-deriv-replace 'mom_yR 'mom_z 'mom_zR)
   'mom_yR
   (flux-deriv-replace 'u 'mom_z 'mom_zR)
   'u
  '(* mom_yR u)
  (symbolic-simp '(- (* mom_yL u) (* mom_yR u)))
   (symbolic-simp-rule '(- (* mom_yL u) (* mom_yR u)))
   '(* (- mom_yL mom_yR) u)
  (symbolic-simp '(* (- mom_yL mom_yR) u))
   (symbolic-simp-rule '(* (- mom_yL mom_yR) u))
   '(* u (- mom_yL mom_yR))
  (symbolic-simp '(* u (- mom_yL mom_yR)))
   (symbolic-simp-rule '(* u (- mom_yL mom_yR)))
    (symbolic-simp-rule 'u)
    'u
    (symbolic-simp-rule '(- mom_yL mom_yR))
     (symbolic-simp-rule 'mom_yL)
     'mom_yL
     (symbolic-simp-rule 'mom_yR)
     'mom_yR
    '(- mom_yL mom_yR)
   '(* u (- mom_yL mom_yR))
  '(* u (- mom_yL mom_yR))
  (flux-deriv-replace '(* mom_z u) 'mom_y 'mom_yL)
   (flux-deriv-replace 'mom_z 'mom_y 'mom_yL)
   'mom_z
   (flux-deriv-replace 'u 'mom_y 'mom_yL)
   'u
  '(* mom_z u)
  (flux-deriv-replace '(* mom_z u) 'mom_z 'mom_zL)
   (flux-deriv-replace 'mom_z 'mom_z 'mom_zL)
   'mom_zL
   (flux-deriv-replace 'u 'mom_z 'mom_zL)
   'u
  '(* mom_zL u)
  (flux-deriv-replace '(* mom_z u) 'mom_y 'mom_yR)
   (flux-deriv-replace 'mom_z 'mom_y 'mom_yR)
   'mom_z
   (flux-deriv-replace 'u 'mom_y 'mom_yR)
   'u
  '(* mom_z u)
  (flux-deriv-replace '(* mom_z u) 'mom_z 'mom_zR)
   (flux-deriv-replace 'mom_z 'mom_z 'mom_zR)
   'mom_zR
   (flux-deriv-replace 'u 'mom_z 'mom_zR)
   'u
  '(* mom_zR u)
  (symbolic-simp '(- (* mom_zL u) (* mom_zR u)))
   (symbolic-simp-rule '(- (* mom_zL u) (* mom_zR u)))
   '(* (- mom_zL mom_zR) u)
  (symbolic-simp '(* (- mom_zL mom_zR) u))
   (symbolic-simp-rule '(* (- mom_zL mom_zR) u))
   '(* u (- mom_zL mom_zR))
  (symbolic-simp '(* u (- mom_zL mom_zR)))
   (symbolic-simp-rule '(* u (- mom_zL mom_zR)))
    (symbolic-simp-rule 'u)
    'u
    (symbolic-simp-rule '(- mom_zL mom_zR))
     (symbolic-simp-rule 'mom_zL)
     'mom_zL
     (symbolic-simp-rule 'mom_zR)
     'mom_zR
    '(- mom_zL mom_zR)
   '(* u (- mom_zL mom_zR))
  '(* u (- mom_zL mom_zR))
  (is-real 0.0 '(mom_y mom_z) '((define u 0.0)))
  #t
  (is-real 0.0 '(mom_y mom_z) '((define u 0.0)))
  #t
  (is-real 0.0 '(mom_y mom_z) '((define u 0.0)))
  #t
 #t
