(function($) {
 Drupal.behaviors.ajax_example = {
 attach:function (context) {
   $('#edit-field-highwire-a-collections-tid').change(function() {
     var tid = $(this).val();
     $.ajax({
       url: '/embo-editors/' + tid,
        success: function(data) {
        $('.editors-info-block').html(data);
       }
     });
   });
 }
}
})(jQuery);
