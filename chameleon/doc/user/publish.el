;; publish.el
;; Emacs publish file for project.
;; Run the following command to execute:
;; emacs --batch --no-init-file --load publish.el --funcall org-publish-all

;; Packages:
(require 'package)
(require 'ox-publish)
(require 'org)
(require 'htmlize)

(setq org-safe-remote-resources '("https://solverstack.gitlabpages.inria.fr/inria-org-html-themes/online-theme-readtheorginria.setup"))

(setq org-html-htmlize-output-type 'css)
(setq org-src-fontify-natively t)
(org-babel-do-load-languages
 'org-babel-load-languages
 '(
   (C . t)
   (fortran . t)
   (python . t)
   (shell . t)
  ))

(setq org-publish-project-alist
      `(("homepage"
         :base-directory "."
         :base-extension "org"
         :publishing-directory "."
         :publishing-function org-html-publish-to-html
         )))
