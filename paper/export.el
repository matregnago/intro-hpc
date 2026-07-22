(add-to-list 'load-path (expand-file-name "."))
(require 'ox-extra)
(ox-extras-activate '(ignore-headlines))
(require 'org-ref)

(add-to-list 'org-latex-packages-alist '("" "url") t)
(add-to-list 'org-latex-packages-alist '("" "sbc-template") t)
(add-to-list 'org-latex-packages-alist '("AUTO" "babel" t ("pdflatex")) t)

(setq enable-local-variables nil
      org-confirm-babel-evaluate nil)

(with-current-buffer (find-file-noselect "sbc.org")
  (org-babel-tangle)
  (org-latex-export-to-latex))
