PRO create_html

  procedures = ['analyze.pro', 'panel.pro', 'rawatom.pro', 'readatom.pro', $
                'readstructure.pro', $
                'viewalpha.pro', 'viewatom.pro', 'viewmolecules.pro', $
                'viewpops.pro', 'viewtermdiag.pro', 'viewtrans.pro']
  mk_html_help, procedures, 'analyze_help.html'
END

