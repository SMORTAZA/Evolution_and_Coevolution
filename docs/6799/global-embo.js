(function ($) {
  //Real Width of Element even if element/parents are hidden
  $.fn.elementRealWidth = function () {
    $clone = this.clone()
        .css("visibility","hidden")
        .appendTo($('body'));
    var $width = $clone.outerWidth();
    $clone.remove();
    return $width;
  };
	
  //Calculate width of text within DOM element
  $.fn.textWidth = function(){
    var eleHtml = $(this).html();
    $(this).html('<span>' + eleHtml + '</span>');
    var eleWidth = $(this).find('span:first').width();
    $(this).html(eleHtml);
    return eleWidth;
  };

	//Check if element is in view area
  $.fn.isIntoView = function () {
    $docViewTop = $(window).scrollTop();
    $docViewBottom = $docViewTop + $(window).height();
    $elemTop = $(this).offset().top;
    $elemBottom = $elemTop + $(this).height();
    $topOffset = $('body').hasClass('admin-menu') ? ($elemTop+100) : ($elemTop+120);
    if($('body').hasClass('not-front')){ 
      $topOffset = $topOffset-10; 
    }
    return (($elemBottom <= $docViewBottom) && ($topOffset >= $docViewTop));
  };
  //Print PDF link using Iframe
  $.fn.printPdf = function (url) {
    var $printWindow = window.open(url, '', 'width=820,height=600');
    $($printWindow).load(function(){
      $printWindow.print();
    });
  }
  //Override aricle author
  Drupal.behaviors.globalFix = {
    attach: function(context, settings) {
      //Fixed floating icons
      $(window).load(function(e) {
        $('.not-front .sidebar-right-wrapper .pane-panels-ajax-tab-tabs li').removeClass('active');

        var menu_wrapper = $('#region-menu');
        var menu_zone = $('#zone-menu');
        var little_menu_offset;

        // SFSF00568055 - now that we have ads, etc in the top we need to
        // calc the menu differently and on page load
        // for pages with a RHS or home carousel, this is an accurate top placement
        if ($('.panel-region-carousel').length) {
          little_menu_offset = $('.panel-region-carousel').offset().top + 10 + 'px';
        }
        else {
          little_menu_offset = $('.sidebar-right-wrapper').offset().top + 'px';
        }
        menu_wrapper.css({'top': little_menu_offset, 'display':'block'});
        
        // onSCroll we need to recalculate, and eventually fix the menu
        $(window).scroll(function(e) {

          if ($('.panel-region-carousel').length) {
            little_menu_offset = $('.panel-region-carousel').offset().top + 10 + 'px';
          }
          else {
            little_menu_offset = $('.sidebar-right-wrapper').offset().top + 'px';
          }
          menu_wrapper.css({'top': little_menu_offset, 'display':'block'});

          if (!menu_zone.isIntoView()) {
            var adminOffset; 
            if ($('body').hasClass('admin-menu')) {
              adminOffset = '35px';
            }
            else {
              adminOffset = '20px'
            }
            menu_wrapper.css({'position':'fixed', 'top':adminOffset});
          }
          else {
            menu_wrapper.css({'top': little_menu_offset, 'display':'block', 'position':'absolute'});
          }
       });
       //Show one accordion by default
       if ($("#highwire_article_accordion_container").length > 0) {
         if($("#highwire_article_accordion_container > .ui-accordion-header").length ==1){
          $( "#highwire_article_accordion_container" ).accordion( "option", "active", 0 ); 
         }          
       }
      });
      $('.not-front .sidebar-right-wrapper .pane-panels-ajax-tab-tabs .panels-ajax-tab-tab').click(function(e){
        $(this).parent('li').addClass('active');
        $('.sidebar-right-wrapper .pane-panels-ajax-tab-tabs + .pane-panels-ajax-tab-container,.sidebar-right-wrapper .pane-panels-ajax-tab-tabs + .panel-separator + .pane-panels-ajax-tab-container').slideDown();
      });
      $('.article-pdf-print').click(function(e){
        e.preventDefault();
        $(this).printPdf($(this).attr('href'));
      });
      $("#region-menu").hover(function () {
        $(this).animate({width: "260px"}, 500 );
        $(this).addClass('hover');
      }, function () {
        $(this).stop(true); //remove the remaining functions/animations
        $(this).removeClass('hover');
        $(this).css({width: "50px"});
        var menu_zone = $('#zone-menu');
        var menu_wrapper = $('#region-menu');
        if(!menu_zone.isIntoView()){
          menu_wrapper.css({'position':'fixed', 'top':'20px'});
          if($('body').hasClass('admin-menu')){
            menu_wrapper.css({'position':'fixed', 'top':'35px'});  
          }
        }
      });
    if($.browser.mozilla || $.browser.msie) {
      $('body').addClass('embo-ie');
      $(".form-type-select .embo-select select").each(function() {
          $width = $(this).width();
          if($width <=0 ){
            $width = $(this).elementRealWidth();
          }else{
            $width = $(this).outerWidth(true);
          }
          $icon_width = $width + 20;
          $(this).css('width',$icon_width);
          $(this).parent('.embo-select').css('width',$width);
          //$(this).siblings('.icon-reorder').css('margin-left', $width);
        });
      }
      if($.browser.mozilla) {
          $('body').addClass('embo-mozilla')
      }
      if ($.browser.msie){//IE fixes
        $('form .button-wrapper .icon-search').click(function(){
          $(this).parents('form:first').submit();
        });
        //Custom Added for Alert Dialog Box
         $(".form-item-frequency .embo-select select").each(function() {
            $width = $(this).outerWidth()+40;
            $icon_width = $width + 20;
            $(this).css('width',$icon_width);
            $(this).parent('.embo-select').css('width',$width);
            $('body').addClass('embo-ie');
        });
      }else{
        if(!$.browser.mozilla)
          $('body').addClass('embo-no-ie');
      }
      
      //Articles page fixes
      $('#fig-data-figures > .fig-data-title-jump .fig-data-title-jump-link:has(h2)').appendTo('#fig-data-figures > .fig-data-title-jump');
      
      //leader ads - repostion if there's no ad
      if ($("#region-ad-top .region-inner a").length == 0) {
        $('body').addClass('top-ads-disabled');
        $('body').css('background-position', '0 0');
        $('.front .page').css('background-position', 'center 85px');
        $('.ad-label').hide();
      } 
      else{
        $('body').addClass('top-ads-enabled');
      }
      $('.highwire-figure.colorbox a.colorbox').removeAttr( "title" );
      $('.highwire-figure.colorbox a.colorbox').colorbox({
        title:function(){
          //$("<div/>").html($(this).data('figure-caption')).text();
          var decodedHTML = $(this).parents('.highwire-figure').next('.fig-caption').html();
          return decodedHTML;
      }}) 
    
      //Article author swap
      if ($('.page-node .highwire-article-citation .highwire-cite-metadata').length > 0) {
        var citeMeta = $('.page-node .highwire-article-citation .highwire-cite-metadata');
        if(! citeMeta.hasClass('loaded')){
          citeMeta.addClass('loaded');
          var citeHtml = citeMeta.html();
          var articleAuthor = citeMeta.parents('.pane-highwire-article-citation').siblings('.author-affiliates');
          articleAuthor.after('<div class="highwire-article-citation highwire-citation-jnl-embo-full-citation"><div class="highwire-cite cite-js"><div class="highwire-cite-metadata"></div></div></div>');
          articleAuthor.next('.highwire-article-citation').find('.highwire-cite-metadata').html(citeHtml);
          citeMeta.html('');
        }        
      }

      //Inline elements
      if ($('.sidebar-right-wrapper .highwire-article-collections > .highwire-list').length > 0) {
        var eleOuterWidth = 230;
        var levelCount = 1;

        //Handle text wrapping with span tag if width exceeds
        $.fn.emboInlineLink = function(level){
          var eleTextWidth = $(this).textWidth();
          if(eleTextWidth >eleOuterWidth){
            $(this).addClass('multiline');
            var eleText = $(this).text();

            if(level){
              var newEleObj = $(this);
            }else{
              $(this).html('<span class="highlight-line line-'+levelCount+'"></span>');
              var newEleObj = $(this).find('.highlight-line.line-'+levelCount);
              newEleObj.html(eleText);
            }

            var newLinetext = new Array();
            var eleTextArr = eleText.trim().split(/\s+/);
            while(eleTextWidth > eleOuterWidth){
              newLinetext.push(eleTextArr.pop());
              newEleObj.html(eleTextArr.join(' '));
              eleTextWidth =  newEleObj.textWidth();
            }
            levelCount = levelCount + 1;
            if(level){
              $(this).after('<span class="highlight-line line-' + levelCount +'"></span>');
              var eleNewLine = $(this).siblings('.line-'+levelCount);
              eleNewLine.html(newLinetext.reverse().join(' '));
            }else{
              $(this).append('<span class="highlight-line line-' + levelCount +'"></span>');
              var eleNewLine = $(this).find('.line-'+levelCount);
              eleNewLine.html(newLinetext.reverse().join(' '));
            }
            eleNewLine.emboInlineLink(true);
          }
          levelCount ++;
        }

        $( ".highwire-article-collections .highwire-list a" ).each(function(){
        $(this).emboInlineLink(false);
        })
      }
    }
  };
})(jQuery);

// Override exposed filter form field's setting on Archive page
jQuery(document).ready(function(){
  jQuery('select#edit-field-highwire-a-collections-tid').val('');
  jQuery('select#edit-field-highwire-article-category-tid').val('');
  jQuery('select#edit-field-highwire-a-epubdate-value-value-year').val('');
  jQuery('select#edit-field-highwire-a-epubdate-value-value-year option:selected').text(Drupal.t('Year'));

  //add new class for empty & CC license on article pages
  (jQuery('body.page-node div.license > p:empty')).parent().addClass('empty');
  (jQuery('body.page-node div.license > p:contains("Creative Commons")')).parent().addClass('by-cc');

  //if a highlight image exists but src contains graphic or mml, the image will be hidden so add class to identify that
  (jQuery('.highwire-cite-highlight img[src*="graphic"]')).parent().parent().addClass('no-display');
  (jQuery('.highwire-cite-highlight img[src*="mml"]')).parent().parent().addClass('no-display');
});
