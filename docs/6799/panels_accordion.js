(function ($) {
  Drupal.behaviors.panels_accordion = {
    attach: function (context, settings) {
      $.each(Drupal.settings.panels_accordion, function($index_id, $options){
        $('#' + $index_id).accordion($options);
      });
    }
  };
}(jQuery));

