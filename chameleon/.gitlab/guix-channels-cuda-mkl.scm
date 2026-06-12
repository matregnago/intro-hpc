;; https://guix.bordeaux.inria.fr/eval/9220640
;; 2026-03-16 10:26:36

(list (channel
        (name 'guix-science-nonfree)
        (url "https://codeberg.org/guix-science/guix-science-nonfree.git")
        (branch "master")
        (commit
          "268e6d26c3e7cd4d0137dc04fe1b0298a4aa2536")
        (introduction
          (make-channel-introduction
            "58661b110325fd5d9b40e6f0177cc486a615817e"
            (openpgp-fingerprint
              "CA4F 8CF4 37D7 478F DA05  5FD4 4213 7701 1A37 8446"))))
      (channel
        (name 'guix)
        (url "https://git.guix.gnu.org/guix.git")
        (branch "master")
        (commit
          "9e17346b663b7b1249c6bb39f1d47a6b0ae7729c")
        (introduction
          (make-channel-introduction
            "9edb3f66fd807b096b48283debdcddccfea34bad"
            (openpgp-fingerprint
              "BBB0 2DDF 2CEA F6A8 0D1D  E643 A2A0 6DF2 A33A 54FA"))))
      (channel
        (name 'guix-science)
        (url "https://codeberg.org/guix-science/guix-science.git")
        (branch "master")
        (commit
          "0cba59f29e1adb47e3b4cca4ec70b0f80fe0e180")
        (introduction
          (make-channel-introduction
            "b1fe5aaff3ab48e798a4cce02f0212bc91f423dc"
            (openpgp-fingerprint
              "CA4F 8CF4 37D7 478F DA05  5FD4 4213 7701 1A37 8446"))))
      (channel
        (name 'guix-past)
        (url "https://codeberg.org/guix-science/guix-past.git")
        (branch "master")
        (commit
          "f5d1cdc58ea810c9d1f1858289b01e9fdc3e4b8f")
        (introduction
          (make-channel-introduction
            "0c119db2ea86a389769f4d2b9c6f5c41c027e336"
            (openpgp-fingerprint
              "3CE4 6455 8A84 FDC6 9DB4  0CFB 090B 1199 3D9A EBB5")))))
