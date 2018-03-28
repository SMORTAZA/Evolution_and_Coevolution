/**
 * @file
 *
 * Behaviors for the article supplementary materials - collapsing and expanding them.
 */

(function () {

  "use strict";

  // Selector for collapsible tags. We can't rely just on .collapsible because
  // it may conflict with Drupal's own collapse functionality.
  var collapsibleTagSelector = '.supplementary-material.collapsible,.type-supplementary-material.collapsible,.section.collapsible,.app.collapsible,.boxed-text.collapsible';

  // Event IDs.
  var eventCollapse = 'collapse';
  var eventExpand = 'expand';

  // Class for the entire toggle
  var wrapClass = 'highwire-collapsible-supplemental';

  // Class for the toggle title.
  var titleClass = 'highwire-collapsible-supplemental-title';
  
  // Class for the toggle content.
  var contentClass = 'highwire-collapsible-supplemental-content';

  // Drupal behavior.
  Drupal.behaviors.HighWire_CollapsibleSupplemental = {
    attach: function (context, settings) {
      applyToggle(context);
    }
  };

  /**
   * Apply toggles to the context. Individuals and global as well.
   *
   * @param context
   *  Drupal behavior context.
   */
  function applyToggle(context) {
    // Settings coming from the backend.
    var titleExpand, titleCollapse, titleAll, expandAll, titleSubtitle;
    if (Drupal.settings.hasOwnProperty('highwire') && Drupal.settings.highwire.hasOwnProperty('collapsible_supplemental')) {
      titleExpand = Drupal.t(Drupal.settings.highwire.collapsible_supplemental.title_expand);
      titleCollapse = Drupal.t(Drupal.settings.highwire.collapsible_supplemental.title_collapse);
      titleAll = Drupal.t(Drupal.settings.highwire.collapsible_supplemental.title_all);
      titleSubtitle = Drupal.settings.highwire.collapsible_supplemental.title_subtitle;
      expandAll = Drupal.settings.highwire.collapsible_supplemental.expand_all;
    }

    // Add individual toggles.
    addToggles(context, titleCollapse, titleExpand, titleSubtitle);
    
    // Main toggle - if needed.
    if (expandAll) {
      addToggleAll(context, titleAll);
    }
  }

  /**
   * Apply individual toggles on the context.
   *
   * @param context
   *  Drupal behavior context.
   * @param titleCollapse
   *  Title for collapsing.
   * @param titleExpand
   *  Title for expanding.
   */
  function addToggles(context, titleCollapse, titleExpand, titleSubtitle) {

    jQuery(collapsibleTagSelector, context).each(function ( idx, item ) {
      var subtitle = '';
      var caption  = '';

      if (titleSubtitle) {
        subtitle = jQuery('.table-label, .fig-label, .collapsible', item).text();
        if (subtitle) {
          subtitle = ' - ' + subtitle;
          caption = jQuery('.caption-title', item).text();
          if (caption) {
            subtitle = subtitle + ' ' + caption;
          }
        }
      }
      var $title = jQuery('<div class="' + titleClass + '">' + titleExpand + subtitle + '</div>');
      var $this = jQuery(item);

      // Attach xref anchor links to collapsible wrapper div and expand on click
      var $xid = $this.attr('id') + '-collapse-anchor';
      jQuery('a[href="#' + $this.attr('id') + '"]').each(function(){
        var link = jQuery(this);
        link.attr('href', '#' + $xid);
        link.bind('click', function(){
          $title.trigger(eventExpand);
        });
      });

      // Wrapper for everything
      $this.wrap('<div id="' + $xid + '" class="'+ wrapClass +'"></div>');
      // Insert Title element
      $this.before($title);
      // Content wrapper
      $this.wrap('<div class="'+ contentClass +'"></div>');

      $title.bind(eventCollapse, function( ) {
        $this.hide();
        $title.html(titleExpand + subtitle);
      });

      $title.bind(eventExpand, function ( ) {
        $this.show();
        $title.html(titleCollapse);
      });

      $title.bind('click', function ( ) {
        $this.is(':visible') ? $title.trigger(eventCollapse) : $title.trigger(eventExpand);
      });

      $title.trigger(eventCollapse);
    });
  }

  /**
   * Creates the toggle-all button and behavior.
   *
   * @param context
   *  Drupal behavior context.
   * @param titleAll
   *  Title for toggling all.
   */
  function addToggleAll(context, titleAll) {
    var $article = jQuery('.article', context);
    var $toggleCount = jQuery(collapsibleTagSelector).length;

    if ($article.length == 0 || $toggleCount == 0) {
      // Return if not the article scope or if there are no collapsible elements
      return;
    }

    var $toggleWrapper = jQuery('<div id="itoggle" class="hw-itoggle-collapsible"><span>' + titleAll + '</span><input type="checkbox" data-label-on="' + Drupal.t('on') + '" data-label-off="' + Drupal.t('off') + '"/></div>')
    $article.prepend($toggleWrapper);

    jQuery('input', $toggleWrapper).iToggle({
      onClickOn: function() {
        jQuery('.' + titleClass).trigger(eventExpand);
      },
      onClickOff: function() {
        jQuery('.' + titleClass).trigger(eventCollapse);
      }
    });
  }

})();
