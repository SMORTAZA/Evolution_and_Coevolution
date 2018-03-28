/**
 * @file
 *
 * Behaviors for generic toggle trigger UI element.
 */
(function ($) {
  Drupal.behaviors.highwireShowNPSSurveyForm = {
    attach: function (context, settings) {
      window.onload = function() {
        var vars = [], arg;
        var args = window.location.href.slice(window.location.href.indexOf('?') + 1).split('&');
        for (var i = 0; i < args.length; i++) {
          arg = args[i].split('=');
          if (arg[0] == 'survey' && arg[1] == '1') {
            (function(t,e,c,n){var o,s,a;t.SMCX=t.SMCX||[],e.getElementById(n)||(o=e.getElementsByTagName(c),s=o[o.length-1],a=e.createElement(c),a.type="text/javascript",a.async=!0,a.id=n,a.src=["https:"===location.protocol?"https://":"http://","widget.surveymonkey.com/collect/website/js/sVMZ3egPEcjKtnfPoVS13yRcF427p8czP1_2BUEi861da7c5ul7EaGM5C_2FGLQkOjar.js"].join(""),s.parentNode.insertBefore(a,s))})(window,document,"script","smcx-sdk");

            return false;
          }
        }
      };
    }
  }
}(jQuery));
