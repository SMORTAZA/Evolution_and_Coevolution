/**
 * Highwire Article Nav
 *
 * Copyright (c) 2010-2011 Board of Trustees, Leland Stanford Jr. University
 * This software is open-source licensed under the GNU Public License Version 2 or later
 * The full license is available in the LICENSE.TXT file at the root of this repository
 */

(function ($) {
  Drupal.behaviors.highwireNavFloat = {
    attach: function (context, settings) {

      /**
       * Function to check if the nav pane should be floated.
       */
      var checkNavFloat = function(threshold, threshold_height, this_pane_height, pane_last_offset, y_scroll_pos) {
        var currentLayout = Drupal.highwireResponsive.getCurrentLayout();
        if ((y_scroll_pos < threshold) || (threshold_height < (this_pane_height + y_scroll_pos)) || (currentLayout == 'mobile')) {
          return false;
        }
        else if (y_scroll_pos > pane_last_offset.top) {
          return true;
        }
      }

      /**
       * Position the floated nav pane.
       */
      $.fn.floatNav = function(options) {
        var settings = $.extend({
          topOffset: 20,
          wrapDiv: '<div class="highwire-nav-float-wrapper"></div>',
        }, options);

        this.wrap(settings.wrapDiv);
        this.parent()
          .hide()
          .css('position', 'fixed')
          .css('top', settings.topOffset + 'px')
          .fadeIn('slow');
      }

      $('.highwire-nav-float', context).once('highwire-nav-float', function() {
        var $wrapper = $(this);
        var $list = $('ul, ol', $wrapper);

        // Scrolling float
        if ($list.data('highwire-float') == '1') {
          var wrap_class = $list.data('highwire-float-class');
          var wrapDiv = '<div class="highwire-nav-float-wrapper ' + wrap_class + '"></div>';
          var topOffset = 20;
          var floated = false;

          // Since we are calculating heights of elements below, should be done after all the elements are loaded
          $(window).load(function(e) {
            var $this_pane = $wrapper.parent().parent();
            var $pane_last = $this_pane.siblings().last();
            var pane_last_offset = $pane_last.offset();
            var pane_last_height = $pane_last.height();
            var threshold = parseInt(pane_last_offset.top) + parseInt(pane_last_height);
            var this_pane_height = $this_pane.height();
            var $footer_section = $('#section-footer');
            var y_scroll_pos = window.pageYOffset;
            var threshold_height = $footer_section.offset().top - topOffset;

            // Check if pane should be floated and set its position.
            check_float_result = checkNavFloat(threshold, threshold_height, this_pane_height, pane_last_offset, y_scroll_pos);
            if (check_float_result === false && floated === true) {
              $this_pane.unwrap();
              floated = false;
            }
            else if (check_float_result === true && floated === false) {
              $this_pane.floatNav({topOffset: topOffset, wrapDiv: wrapDiv});
              floated = true;
            }

            $(window).scroll(function(e) {
              y_scroll_pos = window.pageYOffset;

              pane_last_offset = $pane_last.offset();
              pane_last_height = $pane_last.height();
              threshold = parseInt(pane_last_offset.top) + parseInt(pane_last_height);
              // JCORE-2464:
              // Recalculate threshold height while scrolling due to lazyloaded images.
              threshold_height = $footer_section.offset().top - topOffset;

              // Check if pane should be floated and set its position.
              check_float_result = checkNavFloat(threshold, threshold_height, this_pane_height, pane_last_offset, y_scroll_pos);
              if (check_float_result === false && floated === true) {
                $this_pane.unwrap();
                floated = false;
              }
              else if (check_float_result === true && floated === false) {
                $this_pane.floatNav({topOffset: topOffset, wrapDiv: wrapDiv});
                floated = true;
              }
            });
          });
        }
      });
    }
  };
})(jQuery);
