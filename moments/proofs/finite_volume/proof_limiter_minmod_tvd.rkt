#lang racket

(require "../prover_core.rkt")

 (prove-flux-limiter-tvd '#hash((limiter-expr . (max 0.0 (min 1.0 r))) (limiter-ratio . r) (name . "min_mod")))
  (symbolic-simp '(max 0.0 (min 1.0 r)))
   (symbolic-simp-rule '(max 0.0 (min 1.0 r)))
   '(+ (* 0.5 (+ 0.0 (min 1.0 r))) (* 0.5 (abs (- 0.0 (min 1.0 r)))))
  (symbolic-simp '(+ (* 0.5 (+ 0.0 (min 1.0 r))) (* 0.5 (abs (- 0.0 (min 1.0 r))))))
   (symbolic-simp-rule '(+ (* 0.5 (+ 0.0 (min 1.0 r))) (* 0.5 (abs (- 0.0 (min 1.0 r))))))
    (symbolic-simp-rule '(* 0.5 (+ 0.0 (min 1.0 r))))
    '(+ (* 0.5 0.0) (* 0.5 (min 1.0 r)))
    (symbolic-simp-rule '(* 0.5 (abs (- 0.0 (min 1.0 r)))))
     (symbolic-simp-rule 0.5)
     0.5
     (symbolic-simp-rule '(abs (- 0.0 (min 1.0 r))))
      (symbolic-simp-rule '(- 0.0 (min 1.0 r)))
      '(* -1.0 (min 1.0 r))
     '(abs (* -1.0 (min 1.0 r)))
    '(* 0.5 (abs (* -1.0 (min 1.0 r))))
   '(+ (+ (* 0.5 0.0) (* 0.5 (min 1.0 r))) (* 0.5 (abs (* -1.0 (min 1.0 r)))))
  (symbolic-simp '(+ (+ (* 0.5 0.0) (* 0.5 (min 1.0 r))) (* 0.5 (abs (* -1.0 (min 1.0 r))))))
   (symbolic-simp-rule '(+ (+ (* 0.5 0.0) (* 0.5 (min 1.0 r))) (* 0.5 (abs (* -1.0 (min 1.0 r))))))
   '(+ (* 0.5 0.0) (+ (* 0.5 (min 1.0 r)) (* 0.5 (abs (* -1.0 (min 1.0 r))))))
  (symbolic-simp '(+ (* 0.5 0.0) (+ (* 0.5 (min 1.0 r)) (* 0.5 (abs (* -1.0 (min 1.0 r)))))))
   (symbolic-simp-rule '(+ (* 0.5 0.0) (+ (* 0.5 (min 1.0 r)) (* 0.5 (abs (* -1.0 (min 1.0 r)))))))
    (symbolic-simp-rule '(* 0.5 0.0))
    0.0
    (symbolic-simp-rule '(+ (* 0.5 (min 1.0 r)) (* 0.5 (abs (* -1.0 (min 1.0 r))))))
     (symbolic-simp-rule '(* 0.5 (min 1.0 r)))
      (symbolic-simp-rule 0.5)
      0.5
      (symbolic-simp-rule '(min 1.0 r))
      '(- (* 0.5 (+ 1.0 r)) (* 0.5 (abs (- 1.0 r))))
     '(* 0.5 (- (* 0.5 (+ 1.0 r)) (* 0.5 (abs (- 1.0 r)))))
     (symbolic-simp-rule '(* 0.5 (abs (* -1.0 (min 1.0 r)))))
      (symbolic-simp-rule 0.5)
      0.5
      (symbolic-simp-rule '(abs (* -1.0 (min 1.0 r))))
      '(abs (min 1.0 r))
     '(* 0.5 (abs (min 1.0 r)))
    '(+ (* 0.5 (- (* 0.5 (+ 1.0 r)) (* 0.5 (abs (- 1.0 r))))) (* 0.5 (abs (min 1.0 r))))
   '(+ 0.0 (+ (* 0.5 (- (* 0.5 (+ 1.0 r)) (* 0.5 (abs (- 1.0 r))))) (* 0.5 (abs (min 1.0 r)))))
  (symbolic-simp '(+ 0.0 (+ (* 0.5 (- (* 0.5 (+ 1.0 r)) (* 0.5 (abs (- 1.0 r))))) (* 0.5 (abs (min 1.0 r))))))
   (symbolic-simp-rule '(+ 0.0 (+ (* 0.5 (- (* 0.5 (+ 1.0 r)) (* 0.5 (abs (- 1.0 r))))) (* 0.5 (abs (min 1.0 r))))))
   '(+ (* 0.5 (- (* 0.5 (+ 1.0 r)) (* 0.5 (abs (- 1.0 r))))) (* 0.5 (abs (min 1.0 r))))
  (symbolic-simp '(+ (* 0.5 (- (* 0.5 (+ 1.0 r)) (* 0.5 (abs (- 1.0 r))))) (* 0.5 (abs (min 1.0 r)))))
   (symbolic-simp-rule '(+ (* 0.5 (- (* 0.5 (+ 1.0 r)) (* 0.5 (abs (- 1.0 r))))) (* 0.5 (abs (min 1.0 r)))))
    (symbolic-simp-rule '(* 0.5 (- (* 0.5 (+ 1.0 r)) (* 0.5 (abs (- 1.0 r))))))
     (symbolic-simp-rule 0.5)
     0.5
     (symbolic-simp-rule '(- (* 0.5 (+ 1.0 r)) (* 0.5 (abs (- 1.0 r)))))
     '(* 0.5 (- (+ 1.0 r) (abs (- 1.0 r))))
    '(* 0.5 (* 0.5 (- (+ 1.0 r) (abs (- 1.0 r)))))
    (symbolic-simp-rule '(* 0.5 (abs (min 1.0 r))))
     (symbolic-simp-rule 0.5)
     0.5
     (symbolic-simp-rule '(abs (min 1.0 r)))
      (symbolic-simp-rule '(min 1.0 r))
      '(- (* 0.5 (+ 1.0 r)) (* 0.5 (abs (- 1.0 r))))
     '(abs (- (* 0.5 (+ 1.0 r)) (* 0.5 (abs (- 1.0 r)))))
    '(* 0.5 (abs (- (* 0.5 (+ 1.0 r)) (* 0.5 (abs (- 1.0 r))))))
   '(+ (* 0.5 (* 0.5 (- (+ 1.0 r) (abs (- 1.0 r))))) (* 0.5 (abs (- (* 0.5 (+ 1.0 r)) (* 0.5 (abs (- 1.0 r)))))))
  (symbolic-simp '(+ (* 0.5 (* 0.5 (- (+ 1.0 r) (abs (- 1.0 r))))) (* 0.5 (abs (- (* 0.5 (+ 1.0 r)) (* 0.5 (abs (- 1.0 r))))))))
   (symbolic-simp-rule '(+ (* 0.5 (* 0.5 (- (+ 1.0 r) (abs (- 1.0 r))))) (* 0.5 (abs (- (* 0.5 (+ 1.0 r)) (* 0.5 (abs (- 1.0 r))))))))
    (symbolic-simp-rule '(* 0.5 (* 0.5 (- (+ 1.0 r) (abs (- 1.0 r))))))
    '(* 0.25 (- (+ 1.0 r) (abs (- 1.0 r))))
    (symbolic-simp-rule '(* 0.5 (abs (- (* 0.5 (+ 1.0 r)) (* 0.5 (abs (- 1.0 r)))))))
     (symbolic-simp-rule 0.5)
     0.5
     (symbolic-simp-rule '(abs (- (* 0.5 (+ 1.0 r)) (* 0.5 (abs (- 1.0 r))))))
      (symbolic-simp-rule '(- (* 0.5 (+ 1.0 r)) (* 0.5 (abs (- 1.0 r)))))
      '(* 0.5 (- (+ 1.0 r) (abs (- 1.0 r))))
     '(abs (* 0.5 (- (+ 1.0 r) (abs (- 1.0 r)))))
    '(* 0.5 (abs (* 0.5 (- (+ 1.0 r) (abs (- 1.0 r))))))
   '(+ (* 0.25 (- (+ 1.0 r) (abs (- 1.0 r)))) (* 0.5 (abs (* 0.5 (- (+ 1.0 r) (abs (- 1.0 r)))))))
  (symbolic-simp '(+ (* 0.25 (- (+ 1.0 r) (abs (- 1.0 r)))) (* 0.5 (abs (* 0.5 (- (+ 1.0 r) (abs (- 1.0 r))))))))
   (symbolic-simp-rule '(+ (* 0.25 (- (+ 1.0 r) (abs (- 1.0 r)))) (* 0.5 (abs (* 0.5 (- (+ 1.0 r) (abs (- 1.0 r))))))))
    (symbolic-simp-rule '(* 0.25 (- (+ 1.0 r) (abs (- 1.0 r)))))
     (symbolic-simp-rule 0.25)
     0.25
     (symbolic-simp-rule '(- (+ 1.0 r) (abs (- 1.0 r))))
      (symbolic-simp-rule '(+ 1.0 r))
       (symbolic-simp-rule 1.0)
       1.0
       (symbolic-simp-rule 'r)
       'r
      '(+ 1.0 r)
      (symbolic-simp-rule '(abs (- 1.0 r)))
       (symbolic-simp-rule '(- 1.0 r))
        (symbolic-simp-rule 1.0)
        1.0
        (symbolic-simp-rule 'r)
        'r
       '(- 1.0 r)
      '(abs (- 1.0 r))
     '(- (+ 1.0 r) (abs (- 1.0 r)))
    '(* 0.25 (- (+ 1.0 r) (abs (- 1.0 r))))
    (symbolic-simp-rule '(* 0.5 (abs (* 0.5 (- (+ 1.0 r) (abs (- 1.0 r)))))))
     (symbolic-simp-rule 0.5)
     0.5
     (symbolic-simp-rule '(abs (* 0.5 (- (+ 1.0 r) (abs (- 1.0 r))))))
      (symbolic-simp-rule '(* 0.5 (- (+ 1.0 r) (abs (- 1.0 r)))))
       (symbolic-simp-rule 0.5)
       0.5
       (symbolic-simp-rule '(- (+ 1.0 r) (abs (- 1.0 r))))
        (symbolic-simp-rule '(+ 1.0 r))
         (symbolic-simp-rule 1.0)
         1.0
         (symbolic-simp-rule 'r)
         'r
        '(+ 1.0 r)
        (symbolic-simp-rule '(abs (- 1.0 r)))
         (symbolic-simp-rule '(- 1.0 r))
          (symbolic-simp-rule 1.0)
          1.0
          (symbolic-simp-rule 'r)
          'r
         '(- 1.0 r)
        '(abs (- 1.0 r))
       '(- (+ 1.0 r) (abs (- 1.0 r)))
      '(* 0.5 (- (+ 1.0 r) (abs (- 1.0 r))))
     '(abs (* 0.5 (- (+ 1.0 r) (abs (- 1.0 r)))))
    '(* 0.5 (abs (* 0.5 (- (+ 1.0 r) (abs (- 1.0 r))))))
   '(+ (* 0.25 (- (+ 1.0 r) (abs (- 1.0 r)))) (* 0.5 (abs (* 0.5 (- (+ 1.0 r) (abs (- 1.0 r)))))))
  '(+ (* 0.25 (- (+ 1.0 r) (abs (- 1.0 r)))) (* 0.5 (abs (* 0.5 (- (+ 1.0 r) (abs (- 1.0 r)))))))
  (symbolic-simp-positive '(+ (* 0.25 (- (+ 1.0 r) (abs (- 1.0 r)))) (* 0.5 (abs (* 0.5 (- (+ 1.0 r) (abs (- 1.0 r))))))) 'r)
   (symbolic-simp-positive-rule '(+ (* 0.25 (- (+ 1.0 r) (abs (- 1.0 r)))) (* 0.5 (abs (* 0.5 (- (+ 1.0 r) (abs (- 1.0 r))))))) 'r)
    (symbolic-simp-positive-rule '(* 0.25 (- (+ 1.0 r) (abs (- 1.0 r)))) 'r)
     (symbolic-simp-positive-rule 0.25 'r)
     0.25
     (symbolic-simp-positive-rule '(- (+ 1.0 r) (abs (- 1.0 r))) 'r)
      (symbolic-simp-positive-rule '(+ 1.0 r) 'r)
       (symbolic-simp-positive-rule 1.0 'r)
       1.0
       (symbolic-simp-positive-rule 'r 'r)
       'r
      '(+ 1.0 r)
      (symbolic-simp-positive-rule '(abs (- 1.0 r)) 'r)
      '(- 1.0 r)
     '(- (+ 1.0 r) (- 1.0 r))
    '(* 0.25 (- (+ 1.0 r) (- 1.0 r)))
    (symbolic-simp-positive-rule '(* 0.5 (abs (* 0.5 (- (+ 1.0 r) (abs (- 1.0 r)))))) 'r)
     (symbolic-simp-positive-rule 0.5 'r)
     0.5
     (symbolic-simp-positive-rule '(abs (* 0.5 (- (+ 1.0 r) (abs (- 1.0 r))))) 'r)
     '(* 0.5 (- (+ 1.0 r) (abs (- 1.0 r))))
    '(* 0.5 (* 0.5 (- (+ 1.0 r) (abs (- 1.0 r)))))
   '(+ (* 0.25 (- (+ 1.0 r) (- 1.0 r))) (* 0.5 (* 0.5 (- (+ 1.0 r) (abs (- 1.0 r))))))
  (symbolic-simp-positive '(+ (* 0.25 (- (+ 1.0 r) (- 1.0 r))) (* 0.5 (* 0.5 (- (+ 1.0 r) (abs (- 1.0 r)))))) 'r)
   (symbolic-simp-positive-rule '(+ (* 0.25 (- (+ 1.0 r) (- 1.0 r))) (* 0.5 (* 0.5 (- (+ 1.0 r) (abs (- 1.0 r)))))) 'r)
    (symbolic-simp-positive-rule '(* 0.25 (- (+ 1.0 r) (- 1.0 r))) 'r)
     (symbolic-simp-positive-rule 0.25 'r)
     0.25
     (symbolic-simp-positive-rule '(- (+ 1.0 r) (- 1.0 r)) 'r)
      (symbolic-simp-positive-rule '(+ 1.0 r) 'r)
       (symbolic-simp-positive-rule 1.0 'r)
       1.0
       (symbolic-simp-positive-rule 'r 'r)
       'r
      '(+ 1.0 r)
      (symbolic-simp-positive-rule '(- 1.0 r) 'r)
       (symbolic-simp-positive-rule 1.0 'r)
       1.0
       (symbolic-simp-positive-rule 'r 'r)
       'r
      '(- 1.0 r)
     '(- (+ 1.0 r) (- 1.0 r))
    '(* 0.25 (- (+ 1.0 r) (- 1.0 r)))
    (symbolic-simp-positive-rule '(* 0.5 (* 0.5 (- (+ 1.0 r) (abs (- 1.0 r))))) 'r)
     (symbolic-simp-positive-rule 0.5 'r)
     0.5
     (symbolic-simp-positive-rule '(* 0.5 (- (+ 1.0 r) (abs (- 1.0 r)))) 'r)
      (symbolic-simp-positive-rule 0.5 'r)
      0.5
      (symbolic-simp-positive-rule '(- (+ 1.0 r) (abs (- 1.0 r))) 'r)
       (symbolic-simp-positive-rule '(+ 1.0 r) 'r)
        (symbolic-simp-positive-rule 1.0 'r)
        1.0
        (symbolic-simp-positive-rule 'r 'r)
        'r
       '(+ 1.0 r)
       (symbolic-simp-positive-rule '(abs (- 1.0 r)) 'r)
       '(- 1.0 r)
      '(- (+ 1.0 r) (- 1.0 r))
     '(* 0.5 (- (+ 1.0 r) (- 1.0 r)))
    '(* 0.5 (* 0.5 (- (+ 1.0 r) (- 1.0 r))))
   '(+ (* 0.25 (- (+ 1.0 r) (- 1.0 r))) (* 0.5 (* 0.5 (- (+ 1.0 r) (- 1.0 r)))))
  (symbolic-simp-positive '(+ (* 0.25 (- (+ 1.0 r) (- 1.0 r))) (* 0.5 (* 0.5 (- (+ 1.0 r) (- 1.0 r))))) 'r)
   (symbolic-simp-positive-rule '(+ (* 0.25 (- (+ 1.0 r) (- 1.0 r))) (* 0.5 (* 0.5 (- (+ 1.0 r) (- 1.0 r))))) 'r)
    (symbolic-simp-positive-rule '(* 0.25 (- (+ 1.0 r) (- 1.0 r))) 'r)
     (symbolic-simp-positive-rule 0.25 'r)
     0.25
     (symbolic-simp-positive-rule '(- (+ 1.0 r) (- 1.0 r)) 'r)
      (symbolic-simp-positive-rule '(+ 1.0 r) 'r)
       (symbolic-simp-positive-rule 1.0 'r)
       1.0
       (symbolic-simp-positive-rule 'r 'r)
       'r
      '(+ 1.0 r)
      (symbolic-simp-positive-rule '(- 1.0 r) 'r)
       (symbolic-simp-positive-rule 1.0 'r)
       1.0
       (symbolic-simp-positive-rule 'r 'r)
       'r
      '(- 1.0 r)
     '(- (+ 1.0 r) (- 1.0 r))
    '(* 0.25 (- (+ 1.0 r) (- 1.0 r)))
    (symbolic-simp-positive-rule '(* 0.5 (* 0.5 (- (+ 1.0 r) (- 1.0 r)))) 'r)
     (symbolic-simp-positive-rule 0.5 'r)
     0.5
     (symbolic-simp-positive-rule '(* 0.5 (- (+ 1.0 r) (- 1.0 r))) 'r)
      (symbolic-simp-positive-rule 0.5 'r)
      0.5
      (symbolic-simp-positive-rule '(- (+ 1.0 r) (- 1.0 r)) 'r)
       (symbolic-simp-positive-rule '(+ 1.0 r) 'r)
        (symbolic-simp-positive-rule 1.0 'r)
        1.0
        (symbolic-simp-positive-rule 'r 'r)
        'r
       '(+ 1.0 r)
       (symbolic-simp-positive-rule '(- 1.0 r) 'r)
        (symbolic-simp-positive-rule 1.0 'r)
        1.0
        (symbolic-simp-positive-rule 'r 'r)
        'r
       '(- 1.0 r)
      '(- (+ 1.0 r) (- 1.0 r))
     '(* 0.5 (- (+ 1.0 r) (- 1.0 r)))
    '(* 0.5 (* 0.5 (- (+ 1.0 r) (- 1.0 r))))
   '(+ (* 0.25 (- (+ 1.0 r) (- 1.0 r))) (* 0.5 (* 0.5 (- (+ 1.0 r) (- 1.0 r)))))
  '(+ (* 0.25 (- (+ 1.0 r) (- 1.0 r))) (* 0.5 (* 0.5 (- (+ 1.0 r) (- 1.0 r)))))
  (symbolic-simp '(+ (+ (* 0.0 (- (+ 1.0 r) (- 1.0 r))) (* 0.25 (- (+ 0.0 1.0) (- 0.0 1.0)))) (+ (* 0.0 (* 0.5 (- (+ 1.0 r) (- 1.0 r)))) (* 0.5 (+ (* 0.0 (- (+ 1.0 r) (- 1.0 r))) (* 0.5 (- (+ 0.0 1.0) (- 0.0 1.0))))))))
   (symbolic-simp-rule '(+ (+ (* 0.0 (- (+ 1.0 r) (- 1.0 r))) (* 0.25 (- (+ 0.0 1.0) (- 0.0 1.0)))) (+ (* 0.0 (* 0.5 (- (+ 1.0 r) (- 1.0 r)))) (* 0.5 (+ (* 0.0 (- (+ 1.0 r) (- 1.0 r))) (* 0.5 (- (+ 0.0 1.0) (- 0.0 1.0))))))))
   '(+ (* 0.0 (- (+ 1.0 r) (- 1.0 r))) (+ (* 0.25 (- (+ 0.0 1.0) (- 0.0 1.0))) (+ (* 0.0 (* 0.5 (- (+ 1.0 r) (- 1.0 r)))) (* 0.5 (+ (* 0.0 (- (+ 1.0 r) (- 1.0 r))) (* 0.5 (- (+ 0.0 1.0) (- 0.0 1.0))))))))
  (symbolic-simp '(+ (* 0.0 (- (+ 1.0 r) (- 1.0 r))) (+ (* 0.25 (- (+ 0.0 1.0) (- 0.0 1.0))) (+ (* 0.0 (* 0.5 (- (+ 1.0 r) (- 1.0 r)))) (* 0.5 (+ (* 0.0 (- (+ 1.0 r) (- 1.0 r))) (* 0.5 (- (+ 0.0 1.0) (- 0.0 1.0)))))))))
   (symbolic-simp-rule '(+ (* 0.0 (- (+ 1.0 r) (- 1.0 r))) (+ (* 0.25 (- (+ 0.0 1.0) (- 0.0 1.0))) (+ (* 0.0 (* 0.5 (- (+ 1.0 r) (- 1.0 r)))) (* 0.5 (+ (* 0.0 (- (+ 1.0 r) (- 1.0 r))) (* 0.5 (- (+ 0.0 1.0) (- 0.0 1.0)))))))))
    (symbolic-simp-rule '(* 0.0 (- (+ 1.0 r) (- 1.0 r))))
    0.0
    (symbolic-simp-rule '(+ (* 0.25 (- (+ 0.0 1.0) (- 0.0 1.0))) (+ (* 0.0 (* 0.5 (- (+ 1.0 r) (- 1.0 r)))) (* 0.5 (+ (* 0.0 (- (+ 1.0 r) (- 1.0 r))) (* 0.5 (- (+ 0.0 1.0) (- 0.0 1.0))))))))
     (symbolic-simp-rule '(* 0.25 (- (+ 0.0 1.0) (- 0.0 1.0))))
      (symbolic-simp-rule 0.25)
      0.25
      (symbolic-simp-rule '(- (+ 0.0 1.0) (- 0.0 1.0)))
       (symbolic-simp-rule '(+ 0.0 1.0))
       1.0
       (symbolic-simp-rule '(- 0.0 1.0))
       '(* -1.0 1.0)
      '(- 1.0 (* -1.0 1.0))
     '(* 0.25 (- 1.0 (* -1.0 1.0)))
     (symbolic-simp-rule '(+ (* 0.0 (* 0.5 (- (+ 1.0 r) (- 1.0 r)))) (* 0.5 (+ (* 0.0 (- (+ 1.0 r) (- 1.0 r))) (* 0.5 (- (+ 0.0 1.0) (- 0.0 1.0)))))))
      (symbolic-simp-rule '(* 0.0 (* 0.5 (- (+ 1.0 r) (- 1.0 r)))))
      0.0
      (symbolic-simp-rule '(* 0.5 (+ (* 0.0 (- (+ 1.0 r) (- 1.0 r))) (* 0.5 (- (+ 0.0 1.0) (- 0.0 1.0))))))
      '(+ (* 0.5 (* 0.0 (- (+ 1.0 r) (- 1.0 r)))) (* 0.5 (* 0.5 (- (+ 0.0 1.0) (- 0.0 1.0)))))
     '(+ 0.0 (+ (* 0.5 (* 0.0 (- (+ 1.0 r) (- 1.0 r)))) (* 0.5 (* 0.5 (- (+ 0.0 1.0) (- 0.0 1.0))))))
    '(+ (* 0.25 (- 1.0 (* -1.0 1.0))) (+ 0.0 (+ (* 0.5 (* 0.0 (- (+ 1.0 r) (- 1.0 r)))) (* 0.5 (* 0.5 (- (+ 0.0 1.0) (- 0.0 1.0)))))))
   '(+ 0.0 (+ (* 0.25 (- 1.0 (* -1.0 1.0))) (+ 0.0 (+ (* 0.5 (* 0.0 (- (+ 1.0 r) (- 1.0 r)))) (* 0.5 (* 0.5 (- (+ 0.0 1.0) (- 0.0 1.0))))))))
  (symbolic-simp '(+ 0.0 (+ (* 0.25 (- 1.0 (* -1.0 1.0))) (+ 0.0 (+ (* 0.5 (* 0.0 (- (+ 1.0 r) (- 1.0 r)))) (* 0.5 (* 0.5 (- (+ 0.0 1.0) (- 0.0 1.0)))))))))
   (symbolic-simp-rule '(+ 0.0 (+ (* 0.25 (- 1.0 (* -1.0 1.0))) (+ 0.0 (+ (* 0.5 (* 0.0 (- (+ 1.0 r) (- 1.0 r)))) (* 0.5 (* 0.5 (- (+ 0.0 1.0) (- 0.0 1.0)))))))))
   '(+ (* 0.25 (- 1.0 (* -1.0 1.0))) (+ 0.0 (+ (* 0.5 (* 0.0 (- (+ 1.0 r) (- 1.0 r)))) (* 0.5 (* 0.5 (- (+ 0.0 1.0) (- 0.0 1.0)))))))
  (symbolic-simp '(+ (* 0.25 (- 1.0 (* -1.0 1.0))) (+ 0.0 (+ (* 0.5 (* 0.0 (- (+ 1.0 r) (- 1.0 r)))) (* 0.5 (* 0.5 (- (+ 0.0 1.0) (- 0.0 1.0))))))))
   (symbolic-simp-rule '(+ (* 0.25 (- 1.0 (* -1.0 1.0))) (+ 0.0 (+ (* 0.5 (* 0.0 (- (+ 1.0 r) (- 1.0 r)))) (* 0.5 (* 0.5 (- (+ 0.0 1.0) (- 0.0 1.0))))))))
    (symbolic-simp-rule '(* 0.25 (- 1.0 (* -1.0 1.0))))
     (symbolic-simp-rule 0.25)
     0.25
     (symbolic-simp-rule '(- 1.0 (* -1.0 1.0)))
      (symbolic-simp-rule 1.0)
      1.0
      (symbolic-simp-rule '(* -1.0 1.0))
      -1.0
     '(- 1.0 -1.0)
    '(* 0.25 (- 1.0 -1.0))
    (symbolic-simp-rule '(+ 0.0 (+ (* 0.5 (* 0.0 (- (+ 1.0 r) (- 1.0 r)))) (* 0.5 (* 0.5 (- (+ 0.0 1.0) (- 0.0 1.0)))))))
    '(+ (* 0.5 (* 0.0 (- (+ 1.0 r) (- 1.0 r)))) (* 0.5 (* 0.5 (- (+ 0.0 1.0) (- 0.0 1.0)))))
   '(+ (* 0.25 (- 1.0 -1.0)) (+ (* 0.5 (* 0.0 (- (+ 1.0 r) (- 1.0 r)))) (* 0.5 (* 0.5 (- (+ 0.0 1.0) (- 0.0 1.0))))))
  (symbolic-simp '(+ (* 0.25 (- 1.0 -1.0)) (+ (* 0.5 (* 0.0 (- (+ 1.0 r) (- 1.0 r)))) (* 0.5 (* 0.5 (- (+ 0.0 1.0) (- 0.0 1.0)))))))
   (symbolic-simp-rule '(+ (* 0.25 (- 1.0 -1.0)) (+ (* 0.5 (* 0.0 (- (+ 1.0 r) (- 1.0 r)))) (* 0.5 (* 0.5 (- (+ 0.0 1.0) (- 0.0 1.0)))))))
    (symbolic-simp-rule '(* 0.25 (- 1.0 -1.0)))
     (symbolic-simp-rule 0.25)
     0.25
     (symbolic-simp-rule '(- 1.0 -1.0))
     2.0
    '(* 0.25 2.0)
    (symbolic-simp-rule '(+ (* 0.5 (* 0.0 (- (+ 1.0 r) (- 1.0 r)))) (* 0.5 (* 0.5 (- (+ 0.0 1.0) (- 0.0 1.0))))))
     (symbolic-simp-rule '(* 0.5 (* 0.0 (- (+ 1.0 r) (- 1.0 r)))))
     '(* 0.0 (- (+ 1.0 r) (- 1.0 r)))
     (symbolic-simp-rule '(* 0.5 (* 0.5 (- (+ 0.0 1.0) (- 0.0 1.0)))))
     '(* 0.25 (- (+ 0.0 1.0) (- 0.0 1.0)))
    '(+ (* 0.0 (- (+ 1.0 r) (- 1.0 r))) (* 0.25 (- (+ 0.0 1.0) (- 0.0 1.0))))
   '(+ (* 0.25 2.0) (+ (* 0.0 (- (+ 1.0 r) (- 1.0 r))) (* 0.25 (- (+ 0.0 1.0) (- 0.0 1.0)))))
  (symbolic-simp '(+ (* 0.25 2.0) (+ (* 0.0 (- (+ 1.0 r) (- 1.0 r))) (* 0.25 (- (+ 0.0 1.0) (- 0.0 1.0))))))
   (symbolic-simp-rule '(+ (* 0.25 2.0) (+ (* 0.0 (- (+ 1.0 r) (- 1.0 r))) (* 0.25 (- (+ 0.0 1.0) (- 0.0 1.0))))))
    (symbolic-simp-rule '(* 0.25 2.0))
    0.5
    (symbolic-simp-rule '(+ (* 0.0 (- (+ 1.0 r) (- 1.0 r))) (* 0.25 (- (+ 0.0 1.0) (- 0.0 1.0)))))
     (symbolic-simp-rule '(* 0.0 (- (+ 1.0 r) (- 1.0 r))))
     0.0
     (symbolic-simp-rule '(* 0.25 (- (+ 0.0 1.0) (- 0.0 1.0))))
      (symbolic-simp-rule 0.25)
      0.25
      (symbolic-simp-rule '(- (+ 0.0 1.0) (- 0.0 1.0)))
       (symbolic-simp-rule '(+ 0.0 1.0))
       1.0
       (symbolic-simp-rule '(- 0.0 1.0))
       '(* -1.0 1.0)
      '(- 1.0 (* -1.0 1.0))
     '(* 0.25 (- 1.0 (* -1.0 1.0)))
    '(+ 0.0 (* 0.25 (- 1.0 (* -1.0 1.0))))
   '(+ 0.5 (+ 0.0 (* 0.25 (- 1.0 (* -1.0 1.0)))))
  (symbolic-simp '(+ 0.5 (+ 0.0 (* 0.25 (- 1.0 (* -1.0 1.0))))))
   (symbolic-simp-rule '(+ 0.5 (+ 0.0 (* 0.25 (- 1.0 (* -1.0 1.0))))))
    (symbolic-simp-rule 0.5)
    0.5
    (symbolic-simp-rule '(+ 0.0 (* 0.25 (- 1.0 (* -1.0 1.0)))))
    '(* 0.25 (- 1.0 (* -1.0 1.0)))
   '(+ 0.5 (* 0.25 (- 1.0 (* -1.0 1.0))))
  (symbolic-simp '(+ 0.5 (* 0.25 (- 1.0 (* -1.0 1.0)))))
   (symbolic-simp-rule '(+ 0.5 (* 0.25 (- 1.0 (* -1.0 1.0)))))
    (symbolic-simp-rule 0.5)
    0.5
    (symbolic-simp-rule '(* 0.25 (- 1.0 (* -1.0 1.0))))
     (symbolic-simp-rule 0.25)
     0.25
     (symbolic-simp-rule '(- 1.0 (* -1.0 1.0)))
      (symbolic-simp-rule 1.0)
      1.0
      (symbolic-simp-rule '(* -1.0 1.0))
      -1.0
     '(- 1.0 -1.0)
    '(* 0.25 (- 1.0 -1.0))
   '(+ 0.5 (* 0.25 (- 1.0 -1.0)))
  (symbolic-simp '(+ 0.5 (* 0.25 (- 1.0 -1.0))))
   (symbolic-simp-rule '(+ 0.5 (* 0.25 (- 1.0 -1.0))))
    (symbolic-simp-rule 0.5)
    0.5
    (symbolic-simp-rule '(* 0.25 (- 1.0 -1.0)))
     (symbolic-simp-rule 0.25)
     0.25
     (symbolic-simp-rule '(- 1.0 -1.0))
     2.0
    '(* 0.25 2.0)
   '(+ 0.5 (* 0.25 2.0))
  (symbolic-simp '(+ 0.5 (* 0.25 2.0)))
   (symbolic-simp-rule '(+ 0.5 (* 0.25 2.0)))
    (symbolic-simp-rule 0.5)
    0.5
    (symbolic-simp-rule '(* 0.25 2.0))
    0.5
   '(+ 0.5 0.5)
  (symbolic-simp '(+ 0.5 0.5))
   (symbolic-simp-rule '(+ 0.5 0.5))
   1.0
  (symbolic-simp 1.0)
   (symbolic-simp-rule 1.0)
   1.0
  1.0
  (symbolic-simp-positive 1.0 'r)
   (symbolic-simp-positive-rule 1.0 'r)
   1.0
  1.0
  (symbolic-simp 0.0)
   (symbolic-simp-rule 0.0)
   0.0
  0.0
  (evaluate-limit '(max 0.0 (min 1.0 r)) 'r 1.0)
   (variable-transform '(max 0.0 (min 1.0 r)) 'r 1.0)
    (variable-transform 'max 'r 1.0)
    'max
    (variable-transform 0.0 'r 1.0)
    0.0
    (variable-transform '(min 1.0 r) 'r 1.0)
     (variable-transform 'min 'r 1.0)
     'min
     (variable-transform 1.0 'r 1.0)
     1.0
     (variable-transform 'r 'r 1.0)
     1.0
    '(min 1.0 1.0)
   '(max 0.0 (min 1.0 1.0))
   (evaluate-limit-rule '(max 0.0 (min 1.0 1.0)) 'r 1.0)
    (evaluate-limit-rule 0.0 'r 1.0)
    0.0
    (evaluate-limit-rule '(min 1.0 1.0) 'r 1.0)
    1.0
   '(max 0.0 1.0)
  (evaluate-limit '(max 0.0 1.0) 'r 1.0)
   (variable-transform '(max 0.0 1.0) 'r 1.0)
    (variable-transform 'max 'r 1.0)
    'max
    (variable-transform 0.0 'r 1.0)
    0.0
    (variable-transform 1.0 'r 1.0)
    1.0
   '(max 0.0 1.0)
   (evaluate-limit-rule '(max 0.0 1.0) 'r 1.0)
   1.0
  (evaluate-limit 1.0 'r 1.0)
   (variable-transform 1.0 'r 1.0)
   1.0
   (evaluate-limit-rule 1.0 'r 1.0)
   1.0
  1.0
  (evaluate-limit '(max 0.0 (min 1.0 r)) 'r 0.0)
   (variable-transform '(max 0.0 (min 1.0 r)) 'r 0.0)
    (variable-transform 'max 'r 0.0)
    'max
    (variable-transform 0.0 'r 0.0)
    0.0
    (variable-transform '(min 1.0 r) 'r 0.0)
     (variable-transform 'min 'r 0.0)
     'min
     (variable-transform 1.0 'r 0.0)
     1.0
     (variable-transform 'r 'r 0.0)
     0.0
    '(min 1.0 0.0)
   '(max 0.0 (min 1.0 0.0))
   (evaluate-limit-rule '(max 0.0 (min 1.0 0.0)) 'r 0.0)
    (evaluate-limit-rule 0.0 'r 0.0)
    0.0
    (evaluate-limit-rule '(min 1.0 0.0) 'r 0.0)
    0.0
   '(max 0.0 0.0)
  (evaluate-limit '(max 0.0 0.0) 'r 0.0)
   (variable-transform '(max 0.0 0.0) 'r 0.0)
    (variable-transform 'max 'r 0.0)
    'max
    (variable-transform 0.0 'r 0.0)
    0.0
    (variable-transform 0.0 'r 0.0)
    0.0
   '(max 0.0 0.0)
   (evaluate-limit-rule '(max 0.0 0.0) 'r 0.0)
   0.0
  (evaluate-limit 0.0 'r 0.0)
   (variable-transform 0.0 'r 0.0)
   0.0
   (evaluate-limit-rule 0.0 'r 0.0)
   0.0
  0.0
  (evaluate-limit '(max 0.0 (min 1.0 r)) 'r 2.0)
   (variable-transform '(max 0.0 (min 1.0 r)) 'r 2.0)
    (variable-transform 'max 'r 2.0)
    'max
    (variable-transform 0.0 'r 2.0)
    0.0
    (variable-transform '(min 1.0 r) 'r 2.0)
     (variable-transform 'min 'r 2.0)
     'min
     (variable-transform 1.0 'r 2.0)
     1.0
     (variable-transform 'r 'r 2.0)
     2.0
    '(min 1.0 2.0)
   '(max 0.0 (min 1.0 2.0))
   (evaluate-limit-rule '(max 0.0 (min 1.0 2.0)) 'r 2.0)
    (evaluate-limit-rule 0.0 'r 2.0)
    0.0
    (evaluate-limit-rule '(min 1.0 2.0) 'r 2.0)
    1.0
   '(max 0.0 1.0)
  (evaluate-limit '(max 0.0 1.0) 'r 2.0)
   (variable-transform '(max 0.0 1.0) 'r 2.0)
    (variable-transform 'max 'r 2.0)
    'max
    (variable-transform 0.0 'r 2.0)
    0.0
    (variable-transform 1.0 'r 2.0)
    1.0
   '(max 0.0 1.0)
   (evaluate-limit-rule '(max 0.0 1.0) 'r 2.0)
   1.0
  (evaluate-limit 1.0 'r 2.0)
   (variable-transform 1.0 'r 2.0)
   1.0
   (evaluate-limit-rule 1.0 'r 2.0)
   1.0
  1.0
  (evaluate-limit '(max 0.0 (min 1.0 r)) 'r +inf.0)
   (variable-transform '(max 0.0 (min 1.0 r)) 'r +inf.0)
    (variable-transform 'max 'r +inf.0)
    'max
    (variable-transform 0.0 'r +inf.0)
    0.0
    (variable-transform '(min 1.0 r) 'r +inf.0)
     (variable-transform 'min 'r +inf.0)
     'min
     (variable-transform 1.0 'r +inf.0)
     1.0
     (variable-transform 'r 'r +inf.0)
     +inf.0
    '(min 1.0 +inf.0)
   '(max 0.0 (min 1.0 +inf.0))
   (evaluate-limit-rule '(max 0.0 (min 1.0 +inf.0)) 'r +inf.0)
    (evaluate-limit-rule 0.0 'r +inf.0)
    0.0
    (evaluate-limit-rule '(min 1.0 +inf.0) 'r +inf.0)
    1.0
   '(max 0.0 1.0)
  (evaluate-limit '(max 0.0 1.0) 'r +inf.0)
   (variable-transform '(max 0.0 1.0) 'r +inf.0)
    (variable-transform 'max 'r +inf.0)
    'max
    (variable-transform 0.0 'r +inf.0)
    0.0
    (variable-transform 1.0 'r +inf.0)
    1.0
   '(max 0.0 1.0)
   (evaluate-limit-rule '(max 0.0 1.0) 'r +inf.0)
   1.0
  (evaluate-limit 1.0 'r +inf.0)
   (variable-transform 1.0 'r +inf.0)
   1.0
   (evaluate-limit-rule 1.0 'r +inf.0)
   1.0
  1.0
 #t
