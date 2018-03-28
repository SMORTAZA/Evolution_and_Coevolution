(function ($) {
  Drupal.theme.prototype.isTouchDevice = function () {
    try {
        document.createEvent("TouchEvent");
        return true;
    } catch (e) {
        return false;
    }
  };
  
  /**
   * Provides a always visible scrolbar on touch devices for subject & editor areas.
   */
  Drupal.behaviors.touchFix = {
    attach: function(context, settings) {
      var isTouchDevice = Drupal.theme('isTouchDevice');
      // load nicescroll library only for touch devices - performance
      if(isTouchDevice){ 
        $.getScript( Drupal.settings.jnl_embo_styles.emboStylesJs + "jquery.nicescroll.min.js", function( data, textStatus, jqxhr ) {
          $('.front .panel-region-sidebar-right .pane-meet-editors .snippet-content,.sidebar-right-wrapper .pane-highwire-subject-collections ul.collection').niceScroll({autohidemode:false, 'cursorcolor':'#808080'});
        });
      }
    }
  };
})(jQuery);