/**
 * Highwire AT Symbol
 *
 * Copyright (c) 2010-2011 Board of Trustees, Leland Stanford Jr. University
 * This software is open-source licensed under the GNU Public License Version 2 or later
 * The full license is available in the LICENSE.TXT file at the root of this repository
 */
(function ($) {
  Drupal.behaviors.highwireTablesMarkupProcessor = {
    attach: function(context, settings) {
      $('a.table-expand-inline', context).once('highwireTablesMarkupProcessor', function() {
        $(this, context).each(function() {
          var $caption, captionHTML, tableHTML;
          var self = this;
          var is_collapsed = true;

          // find the tables those does not have captions :)
          $(self, context).each(function(){
            var $caption_temp;
            $caption_temp = $(self).closest('.table');

            if($caption_temp.find('.table-caption').length == 0) {
              $caption_temp.append( "<div class='table-caption'> &nbsp;</div>" );
              Drupal.attachBehaviors($caption_temp[0]);
            }
          });

          $(self).click(function(event) {
            event.preventDefault();
            if (is_collapsed) {
              // Store the table element and caption.
              $caption = $(self).closest('.table').find('.table-caption');
              captionHTML = $caption.html();

              // Load and store the table markup if none exists.
              if (!tableHTML) {
                $.get($(this).data('table-url'), function(data) {
                  $caption.html(tableHTML = data);
                })
                .fail(function() {
                  $caption.html(captionHTML + '<strong>Sorry, there was a problem fetching the table. Please try again later.</strong>');
                })
                .always(function() {
                  Drupal.attachBehaviors($caption[0]);
                  // Update UI after ajax call - DRQUEST-1335.
                  $(self).closest('.table').toggleClass('table-expand-inline');
                  $(self).html('Collapse inline');
                  is_collapsed = false;
                });
              }

              $(self).html('Collapse inline');
              $(self).closest('.table').find('table').wrap('<div class="table-wrapper"></div>');  // #JCORE-1937

              $caption.html(tableHTML);
              Drupal.attachBehaviors($caption[0]);

              var alt = $(self).closest('.table').find('.table-label').text();
              $(self).closest('.table').find('table').attr('alt', alt);
            }
            else {
              $caption.html(captionHTML);
              $(self).html('View inline');
            }

            $(self).closest('.table').toggleClass('table-expand-inline');
            is_collapsed = !is_collapsed;
          });
        });
        /**
         * Added this colorbox calling function as AJAX tabs naigation holds
         * poping up data into model
         */
        $('a.table-expand-popup', context).each(function() {
          cbsettings = $.extend(settings.colorbox, {title: false});
          $(this).colorbox(cbsettings);
        });
      });
    }
  };

  // Attach drupal behaviors to colorbox loading
  $(document).bind('cbox_complete', function() {
    if (colorbox) {
      Drupal.attachBehaviors($(colorbox)[0]);
    }
  });
})(jQuery);
