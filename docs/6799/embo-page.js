(function ($) {
  /**
   * Provides a always visible scrolbar on touch devices for subject & editor areas.
   */
  Drupal.behaviors.emboPage = {
    attach: function (context, settings) {
      $('.embo-toggle-title + .embo-toggle-desc:not(.opened)')
              .addClass('closed')
              .slideUp()
              .prev('.embo-toggle-title')
              .addClass('closed')
              .removeClass('opened');
      $('.embo-toggle-title + .embo-toggle-desc.opened')
              .prev('.embo-toggle-title')
              .addClass('opened')
              .removeClass('closed');
      $('.embo-toggle-title').click(function () {
        if ($(this).next('.embo-toggle-desc').hasClass('closed')) {
          $(this).next('.embo-toggle-desc').slideDown('slow');
          $(this).next('.embo-toggle-desc').addClass('opened').removeClass('closed');
          $(this).addClass('opened').removeClass('closed');
        } else {
          $(this).next('.embo-toggle-desc').slideUp('slow');
          $(this).next('.embo-toggle-desc').addClass('closed').removeClass('opened');
          $(this).addClass('closed').removeClass('opened');
        }
      });
      $( ".embo-toggle-title a" ).click(function( event ) {
        event.preventDefault();
      });
    }
  };
})(jQuery);