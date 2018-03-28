/**
 * HW Jaavascript Content
 *
 * Copyright (c) 2010-2011 Board of Trustees, Leland Stanford Jr. University
 * This software is open-source licensed under the GNU Public License Version 2 or later
 * The full license is available in the LICENSE.TXT file at the root of this repository
 */

//NOTES - these selectors had to be changed for the js to work with the drupal markup!
//-----------------------------------------------------------------------------------
// old-selector		->	new-selector
// #content-block -> 	#panel-content-block
// #slugline 			-> 	.hw-article-citation

/* jQuery javascript functions for content page */
var gIsFrameset = false;
var gIsOMContent;
var gAuthorAffilList = "div.article div.contributors p.affiliation-list-reveal a.view-more";
var gAuthorNotes = "div.article div.contributors p.author-notes-reveal a.view-more";
var gAuthorContributorList = "div.article div.contributors ol.contributor-list a[href^='#aff']";

var gColMain, gCol2, gCol3, gIsMultiCol, gColTempResize;
var gColOrigHeights = { main: 0, col2: 0, col3: 0 };
var gColLows = { main: 0, col2: 0, col3: 0 };
var isIE6 = (jQuery.browser.msie && (jQuery.browser.version.substr(0,2) == '6.'));

/*['abstract','extract','excerpt']*/
var gSiteOptions = {
	popupViews: ["abstract"],
	popabsParam: "abspop=1",
	refPopupLinkTypes: ["medline"],
	refAbsPopupTimeout: 400,
	expandString: "+",
	contractString: "-",
	suppressDockedNav: true,
	suppressDockedSearchNav: true,
	collapsibleLabels: ["h3","h4","span"],
	isiLinkString: "Loading Web of Science citing article data...",
    'bam-ads': {
        'custom-css':
            {
            //'leaderboard': '/publisher/css/hw-publisher-leaderboard.css'
            }
    }
};
var gSiteOptionUnknown = 'UNKNOWN';

(function ($) {
	Drupal.behaviors.highwire_content = {
  	attach: function(context, settings) {

			if (!(getSiteOption("noAuthAffilCollapse") == true)) {
				handleAuthAffil(gAuthorAffilList,gAuthorNotes,gAuthorContributorList);
			}

			fixWackyReflinksMarkup();

			var defaultDockedNavRules = [
				'', '$(.pane-highwire-article-views)',
				'', '$(.pane-highwire-article-citation)',
				'', '$(.pane-highwire-article-nav)'
			];

			//This conflicts with colorbox...do we really need it?
			if (!(getSiteOption("suppressDockedNav") == true)) {
				setupDockBlock(2, 'docked-nav', 'dockblock', defaultDockedNavRules);
			}



		$('#panel-content-block').prepend('<div id="print-slug" class="print-only"></div>');
		var loc=location.hostname;
		if (document.getElementsByName) {
		  var metaArray = document.getElementsByName('citation_journal_title');
		  for (var i=0; i<metaArray.length; i++) {
		    $('<span class="jnl-title">'+ metaArray[i].content + '</span><span class="jnl-url">'+ loc + '</span>' ).appendTo('#print-slug');
		  }
		};


		$('div#col-2  div.cb-slug').clone().appendTo('#print-slug');


		if ((typeof(gIsFrameset) == 'undefined') || (!gIsFrameset)) {
		if(top!= self) { top.location.href = self.location.href };
	}




	/* make links with rel="external-nw" target a new window */
	$("a[rel*=external-nw]").each(
		function(i) {
			linkNewWin($(this), this);
		}
	);

	/* attach onFocus handlers to quick-search fields, to clear the
	   field when it receives focus IF it contains the default
	   value (such as "Search the site") */
	updateFormInput("#header-qs-search-label", "#header-qs-input", '', '');

	/* handle collapsing content-box areas, with user prefs */
	setupCollapsibles();

	/* handle "most" boxes */
	handleMostBoxes();



	}};
	})(jQuery);


function handleAuthAffil(authorAffilList,authorNotes,authorContributorList) {
	var authAffilMatch = getSiteOption('authAffilMatch','div.article div.contributors ol.affiliation-list:has(li)');
	var authAffil = (authAffilMatch != undefined) ? jQuery(authAffilMatch) : '';
	var disableIfMultAffils = getSiteOption('authAffilDisableMultipleMatches',false);
	if (authAffil.length && ((authAffil.length <= 1) || (!disableIfMultAffils))) {
		var expandStr = getSiteOption('authExpandString', null);
		if (expandStr == null) {
			expandStr = getSiteOption('expandString', '+');
		}
		var contractStr = getSiteOption('authContractString', null);
		if (contractStr == null) {
			contractStr = getSiteOption('contractString', '-');
		}
		// var newP = '<p class="affiliation-list-reveal"><a href="#" class="view-more">' + expandStr + '</a> Author Affiliations</p>';
		// /* add auth affil show/hide p */
		// var contribLists = jQuery( "div.article div.contributors:last ol.contributor-list:has(li)");
		// if (contribLists.length) {
		// 	contribLists.after(newP);
		// }
		// /* hide author affiliations until requested */
		// if (authAffil.length) {
		// 	authAffil.each(
		// 		function (i) {
		// 			modClass(authAffil.eq(i),'hideaffil','showaffil');
		// 		}
		// 	);
		// }
		// jQuery(authorAffilList).click(
		// 	function(e) {
		// 			AuthClickOnAffilButton(authorAffilList,authAffilMatch,expandStr,contractStr,e);
		// 	}
		// );
		/* show author affiliations when affil link is selected */
		jQuery(authorContributorList).click(
			function(e) {
					AuthClickOnContribAffilLink(authorAffilList,authAffilMatch,expandStr,contractStr,e);
			}
		);

	var authNotesMatch = getSiteOption('authNotesMatch','div.article div.contributors ul.author-notes:has(li)');
	var authNotes = (authNotesMatch != undefined) ? jQuery(authNotesMatch) : '';
	var disableIfMultNotes = getSiteOption('authNotesDisableMultipleMatches',false);
	if (authNotes.length && ((authNotes.length <= 1) || (!disableIfMultNotes))) {
		var expandStr = getSiteOption('authExpandString', null);
		if (expandStr == null) {
			expandStr = getSiteOption('expandString', '+');
		}
		var contractStr = getSiteOption('authContractString', null);
		if (contractStr == null) {
			contractStr = getSiteOption('contractString', '-');
		}
        var authAffilList = jQuery("div.article div.contributors ol.affiliation-list:has(li):last"); //get the affiliate-list
        var newNotesP = '<p class="author-notes-reveal"><a href="#" class="view-more">' + expandStr + '</a> Author Notes</p>';
		/* add auth notes show/hide p, and move author-notes after it */
	    authAffilList.addClass("has-authnotes").after(authNotes).after(newNotesP);
		/* hide author notes until requested */
		authNotes.each( //already checked if authNotes present
			function (i) {
				modClass(authNotes.eq(i),'hidenotes','shownotes');
			}
		);
		jQuery(authorNotes).click(
			function(e) {
					AuthClickOnNotesButton(authorNotes,authNotesMatch,expandStr,contractStr,e);
			}
		);
    }
        fixColHeights(1);
	}
}

function figExpandInline() {
	var figlinks = jQuery("div.fig-inline:not(.video-inline) a[href*='expansion']");
	if (figlinks.length) {
		figlinks.each(
			function() {
				var $this = jQuery(this);
				var classAttr = $this.attr("class");
				if (!(classAttr && ((classAttr == 'in-nw') || (classAttr == 'ppt-landing')))) {

                    $this.addClass("fig-inline-link");

					if ($this.text().indexOf('n this window') >= 0) {
						$this.text("In this page");
					}
					var parentDiv = $this.parents("div.fig-inline");
					var href = $this.attr("href");
					jQuery( this).click(
						function(e) {
							swapFig(e, href, parentDiv);
                            parentDiv.find('a.fig-inline-link, a.in-nw-vis').unbind('click');
						}
					);
				}
			}
		);
	}
}




function swapFig(e, href, figWrapperEl) {
	var host = document.location.protocol + "//" + document.location.host;
	var path = document.location.pathname;
	var pathseg = path.substring(0, path.lastIndexOf('/'));
	//var baseAjaxUrl = host + pathseg + '/' + href;
	var baseAjaxUrl;
	if (href.indexOf('http:' == 0)) {
		baseAjaxUrl = href;
	}
	else {
		baseAjaxUrl = host + pathseg + '/' + href;
	}
	var ajaxUrl = baseAjaxUrl + ((href.indexOf('?') >= 0) ? '&' : '?') + 'baseURI=' + ((baseAjaxUrl.indexOf('?') > 0) ? baseAjaxUrl.substring(0, baseAjaxUrl.indexOf('?')): baseAjaxUrl);
	//var ajaxUrl = baseAjaxUrl + ((href.indexOf('?') >= 0) ? '&' : '?') + 'baseURI=' + baseAjaxUrl;
	$.ajax({
		url: ajaxUrl,
		dataType: "html",
		type: "GET",
		error: ajaxErr,
		beforeSend: addFigHeaders,
		success: function(xhtml) {
			addFig(xhtml, figWrapperEl);
		},
		complete: ajaxComplete
	});


	e.preventDefault();
}
function addFigHeaders(req) {
	addCommonHeaders(req);
	addPartHeaders(req);
}

function addFig(xhtmlData, figWrapperEl) {
	if (xhtmlData && !(xhtmlData.indexOf('<html') >= 0)) {
        figWrapperEl.addClass("inline-expansion");

        // pick the replacement image out of the div we get back - there should be only one
        largerImage = jQuery( "img", xhtmlData).filter(":first");
        largerImage.addClass("replaced-figure");

        // get the current image, mark it and hide it
        previousImage = jQuery( "img", figWrapperEl).filter(":first");
        previousImage.addClass("previous-figure");
        previousImage.hide();

        // swap previous image out for larger image
        largerImage.appendTo(previousImage.parent());

        // Remove link to "display in this window" by looking for the link and hiding it's parent
        figWrapperEl.find(".callout .callout-links a.fig-inline-link").parent('li').hide();

		newWindowTargets();
		var lookupRule = "div.article div.fig img";
		var numImagesToLoad = checkUnloadedImgs(lookupRule);
		setTimeout("fixHeightForImages(1" + "," + numImagesToLoad + ",'" + lookupRule + "')", 1000);
		fixColHeights(1);
	}
}

function refLinksNewWindowTarget() {
	jQuery( "div.ref-list div.cit-extra a").each(
		function(i) {
			var origTitle = jQuery(this).attr("title");
			var newTitle = '';
			if ((origTitle == undefined) || (!origTitle)) {
				origTitle = '';
			}
			else {
				newTitle = origTitle + ' ';
			}
			newTitle += '[opens in a new window]';
			jQuery( this).attr("target", "_blank").attr("title", newTitle);
		}
	);
}

function addRefPops() {
	var numMissed = 0;
	var maxToSkip = getSiteOption('refPopsMaxToSkip', 10);
	var i = 1;
	var idroot = "#xref-ref-";
	var el = jQuery( idroot + i + '-1');
	/* not all refs appear in text; compensate */
	while (numMissed < maxToSkip) {
		if (el.length) {
			numMissed = 0;
		el.hover(dispRef, hideRef);
		if ((getSiteOption("isIosRefPops") == true)) {
 			//Also set the ios version of the function
			el.hover(iosdispRef, hideRef);
		}
		var j = 2;
		var el2 = jQuery( idroot + i + "-" + j);
		while (el2.length) {
			el2.hover(dispRef, hideRef);
			if ((getSiteOption("isIosRefPops") == true)) {
 				//Also set the ios version of the function
				el2.hover(iosdispRef, hideRef);
			}
			j++;
			el2 = jQuery( idroot + i + "-" + j);
		}
		}
		else {
			numMissed++;
		}
		i++;
		el = jQuery( idroot + i + "-1");
	}
}

function dispRef(e) {
	var link = jQuery( this).attr("href");
	if(jQuery( "div#hovering-ref").length) {
		jQuery( "div#hovering-ref").remove();
		//alert("hovering-ref div removed on new hover!");
	}
	var linkEl = jQuery( link);
	if (linkEl.length) {
		var citHtml = linkEl.next("div").children("div.cit-metadata");
		if (!(citHtml.length)) {
			citHtml = linkEl.parent().next("div").children("div.cit-metadata");
		}
		if (citHtml.length) {
			var newDiv = '<div id="hovering-ref">' + (citHtml.clone().html()) + '</div>';
			jQuery( "body").append(newDiv);
			var elH = getObjHeight(jQuery( "div#hovering-ref"));
			if ((getSiteOption("isIosRefPops") == true)) {
				jQuery( "div#hovering-ref").css("left", 10).css("top", e.pageY-elH).css("position", "absolute");
			}
			else {
				jQuery( "div#hovering-ref").css("left", e.pageX+10).css("top", e.pageY-elH).css("position", "absolute");
			}
		}
	}
}

function hideRef(e) {
	if(jQuery( "div#hovering-ref").length) {
		jQuery( "div#hovering-ref").remove();
	}
}

function getHWCiting() {
	var citingA = jQuery( "#cb-hw-citing-articles");
	if (citingA.length) {
		var newA = '<a id="cb-loading-hw-cited" href="#">Loading citing article data...</a>';
		citingA.replaceWith(newA);
		var href = citingA.attr("href");
		var id = '';
		if (href && (href.indexOf('?') > 0)) {
			var args = href.substring(href.indexOf('?') + 1).split('&');
			for (var i = 0; i < args.length; i++) {
				if (args[i].toLowerCase().indexOf('legid=') == 0) {
					id = args[i].substring(args[i].indexOf('=') + 1);
					if (id.indexOf('#') > 0) {
						id = id.substring(0, id.indexOf('#'));
					}
				}
			}
			if (!(id == '')) {
				var host = document.location.protocol + "//" + document.location.host;
				var ajaxUrl = host + '/cited-by/' + id.replace(/;/,'/');
				$.ajax({
					url: ajaxUrl,
					dataType: "html",
					type: "GET",
					error: ajaxErr,
					success: addHWCiting,
					complete: ajaxCompleteCitedBy
				});

			}
		}
	}
}
function getHWRelatedURLs() {
	var relatedURLsA = jQuery( "#cb-related-urls");
	var relatedURLsMsg = getSiteOption('relatedWebPageLoadingText','Loading related web page data...');
	if (relatedURLsA.length) {
		var newA = '<a id="cb-loading-related-urls" href="#">' + relatedURLsMsg + '</a>';
		relatedURLsA.replaceWith(newA);
		var href = relatedURLsA.attr("href");
		var id = '';
		if (href && (href.indexOf('?') > 0)) {
			var args = href.substring(href.indexOf('?') + 1).split('&');
			for (var i = 0; i < args.length; i++) {
				if (args[i].toLowerCase().indexOf('legid=') == 0) {
					id = args[i].substring(args[i].indexOf('=') + 1);
					if (id.indexOf('#') > 0) {
						id = id.substring(0, id.indexOf('#'));
					}
				}
			}
			if (!(id == '')) {
				var host = document.location.protocol + "//" + document.location.host;
				var ajaxUrl = host + '/related-web-pages/' + id.replace(/;/,'/');
				$.ajax({
					url: ajaxUrl,
					dataType: "html",
					type: "GET",
					error: ajaxErr,
					success: addRelatedURLs,
					complete: ajaxComplete
				});

			}
		}
	}
}
function getPatientInformData() {
	var pInformA = jQuery( "#cb-patientinform");
	if (pInformA.length) {
		var newA = '<a id="cb-loading-patientinform" href="#">Loading <em>patient</em>INFORMation...</a>';
		pInformA.replaceWith(newA);
		var href = pInformA.attr("href");
		var id = '';
		if (href && (href.indexOf('?') > 0)) {
			var args = href.substring(href.indexOf('?') + 1).split('&');
			for (var i = 0; i < args.length; i++) {
				if (args[i].toLowerCase().indexOf('legid=') == 0) {
					id = args[i].substring(args[i].indexOf('=') + 1);
					if (id.indexOf('#') > 0) {
						id = id.substring(0, id.indexOf('#'));
					}
				}
			}
			if (!(id == '')) {
				var host = document.location.protocol + "//" + document.location.host;
				var ajaxUrl = host + '/related-web-pages/patientinform/' + id.replace(/;/,'/');
				$.ajax({
					url: ajaxUrl,
					dataType: "html",
					type: "GET",
					error: ajaxErr,
					success: addPatientInform,
					complete: ajaxComplete
				});

			}
		}
	}
}
function getISIRelated() {
	var relatedA = jQuery( "#cb-isi-similar-articles");
	if (relatedA.length) {
		var newA = '<a id="cb-isi-similar-articles" href="#">Loading Web of Science article data...</a>';
		relatedA.replaceWith(newA);
		var href = relatedA.attr("href");
		var id = '';
		if (href) {
			var hrefDec = decodeURI(href);
			if ((hrefDec.indexOf('?') > 0)) {
				var args = hrefDec.substring(hrefDec.indexOf('?') + 1).split('&');
				for (var i = 0; i < args.length; i++) {
					var argDec = decodeURIComponent(args[i]);
					if (argDec.toLowerCase().indexOf('access_num=') == 0) {
						id = argDec.substring(argDec.indexOf('=') + 1);
						if (id.indexOf('#') > 0) {
							id = id.substring(0, id.indexOf('#'));
						}
					}
				}
				if (!(id == '')) {
					var host = document.location.protocol + "//" + document.location.host;
					var ajaxUrl = host + '/isi-links/has-related/' + id.replace(/;/,'/');
					$.ajax({
						url: ajaxUrl,
						dataType: "html",
						type: "GET",
						error: ajaxErr,
						success: addISIRelated,
						complete: ajaxComplete
					});
				}
			}
		}
	}
}
function getCiting(service, msg, pathseg, successFn, addlParamString) {
	if (typeof(addlParamString) == "undefined") {
		addlParamString = '';
	}
	var citingA = jQuery( "#cb-" + service + "-citing-articles");
	if (citingA.length) {
		var newA = '<a id="cb-loading-' + service + '-cited" href="#">' + msg + '</a>';
		citingA.replaceWith(newA);
		var href = citingA.attr("href");
		var id = '';
		if (href) {
			var hrefDec = decodeURI(href);
			if ((hrefDec.indexOf('?') > 0)) {
				var args = hrefDec.substring(hrefDec.indexOf('?') + 1).split('&');
				for (var i = 0; i < args.length; i++) {
					var argDec = decodeURIComponent(args[i]);
					if (argDec.toLowerCase().indexOf('access_num=') == 0) {
						id = argDec.substring(argDec.indexOf('=') + 1);
						if (id.indexOf('#') > 0) {
							id = id.substring(0, id.indexOf('#'));
						}
					}
				}
				if (!(id == '')) {
					var host = document.location.protocol + "//" + document.location.host;
					var ajaxUrl = host + '/' + service + '-links/' + pathseg + id.replace(/;/,'/') + (addlParamString != '' ? '?' + addlParamString : '');
					$.ajax({
						url: ajaxUrl,
						dataType: "html",
						type: "GET",
						error: ajaxErr,
						success: successFn,
						complete: ajaxComplete
					});
				}
			}
		}
	}
}

function getISICiting() {
	var citingA = jQuery( "#cb-isi-citing-articles");
	if (citingA.length) {
		var newA = '<a id="cb-loading-isi-cited" href="#">Loading Web of Science citing article data...</a>';
		citingA.replaceWith(newA);
		var href = citingA.attr("href");
		var id = '';
		if (href) {
			var hrefDec = decodeURI(href);
			if ((hrefDec.indexOf('?') > 0)) {
				var args = hrefDec.substring(hrefDec.indexOf('?') + 1).split('&');
				for (var i = 0; i < args.length; i++) {
					var argDec = decodeURIComponent(args[i]);
					if (argDec.toLowerCase().indexOf('access_num=') == 0) {
						id = argDec.substring(argDec.indexOf('=') + 1);
						if (id.indexOf('#') > 0) {
							id = id.substring(0, id.indexOf('#'));
						}
					}
				}
				if (!(id == '')) {
					var host = document.location.protocol + "//" + document.location.host;
					var ajaxUrl = host + '/isi-links/' + id.replace(/;/,'/');
					$.ajax({
						url: ajaxUrl,
						dataType: "html",
						type: "GET",
						error: ajaxErr,
						success: addISICiting,
						complete: ajaxComplete
					});
				}
			}
		}
	}
}
function getEntrezLinks() {
	var entrezDiv = jQuery( "#cb-entrez-links-placeholder");
	if (entrezDiv.length) {
		var entrezA = entrezDiv.children("a");
		if (entrezA) {
			var host = document.location.protocol + "//" + document.location.host;
			var ajaxUrl = host + entrezA.attr("href");
			$.ajax({
				url: ajaxUrl,
				dataType: "html",
				type: "GET",
				error: ajaxErr,
				success: addEntrezLinks,
				complete: ajaxComplete
			});
		}
	}
}

function ajaxErr(req, msg, e) {
}
function ajaxComplete(req, msg) {
}
function ajaxCompleteCitedBy(req, msg) {
}

function updateCBItem(cbItem, newHTML, hasData) {
	var parentItem = cbItem.parents("li").eq(0);
	cbItem.replaceWith(newHTML);
	if (!hasData) {
		// hide the parent li
		if (parentItem.length) {
			modClass(parentItem,"nodata","");
			// check if there are any siblings still being displayed
			var otherItems = parentItem.siblings();
			var allItemsEmpty;
			if (otherItems.length) {
				if (otherItems.length == otherItems.filter(".nodata").length) {
					allItemsEmpty = true
				}
				else {
					allItemsEmpty = false;
				}
			}
			else {
				allItemsEmpty = true;
			}
			if (allItemsEmpty) {
				var cbsection = parentItem.parents("div.cb-section").eq(0);
				if (cbsection.length) {
					modClass(cbsection,"nodata","");
				}
				// do we need to look further?
				if (parentItem.parents("div.cb-section").length > 1) {
					var cbSectionSibs =  cbsection.siblings("div.cb-section");
					if (cbSectionSibs.length) {
						if (cbSectionSibs.length == cbSectionSibs.filter(".nodata").length) {
							allItemsEmpty = true
						}
						else {
							allItemsEmpty = false;
						}
					}
					else {
						allItemsEmpty = true;
					}
					if (allItemsEmpty) {
						var cbgrandsection = parentItem.parents("div.cb-section").eq(1);
						if (cbgrandsection.length) {
							modClass(cbgrandsection,"nodata","");
						}
					}
				}
			}
		}
	}
	// in frameset fix targets on child links, forms
	fixFrameLinks(parentItem.find("a,form"));
	fixColHeights(2);
}
function fixFrameLinks(jqItems) {
	if ((gIsFrameset != null) && gIsFrameset) {
		jqItems.each(
			function(i) {
				var href = jQuery( this).attr("href");
				var action = jQuery( this).attr("action"); // if form
				if ((href != null) || (action != null)) {
					var inFrameAnchor = ((href != null) && (((/frameset=/.test(href)) && (/#/.test(href))) || (href.substring(0,1) == '#')));
					if ((!inFrameAnchor) || (action != null)) {
						if ((navigator.userAgent.indexOf("Firefox") >= 0) && (jQuery( this).hasClass("pdf-direct-link"))) {
							jQuery( this).attr("target", "_blank");
						} else if (getSiteOption("hasFrameLinkTargetFunction", false)) {
							jQuery( this).attr("target", setFrameLinkTarget(jQuery( this)));
						} else {
							jQuery( this).attr("target", "_top");
						}
					}
				}
			}
		);
	}
}
function addRelatedURLs(xhtmlData) {
	var cbA = jQuery( "#cb-loading-related-urls");
	if (gIsFrameset) {
		if (cbA.length) {
			updateCBItem(cbA, '<div id="cb-loaded-related-urls-none">Not available in this view</div>', false);
		}
	}
	else if (xhtmlData && !(xhtmlData.indexOf('<span') >= 0)) {
		jQuery( "#related-urls").replaceWith(xhtmlData);
		var relatedWebPagesLabel = getSiteOption('relatedWebPagesLabel', 'Related Web Pages');
		fixColHeights(1);
		if (cbA.length) {
			updateCBItem(cbA, '<a href="#related-urls">' + relatedWebPagesLabel + '</a>', true);
		}
	}
	else {
		if (cbA.length) {
			updateCBItem(cbA, '<div id="cb-loaded-related-urls-none">No related web pages</div>', false);
		}
	}
}
function addPatientInform(xhtmlData) {
	var cbA = jQuery( "#cb-loading-patientinform");
	if (gIsFrameset) {
		if (cbA.length) {
			updateCBItem(cbA, '<div id="cb-loaded-patientinform-none">Not available in this view</div>', false);
		}
	}
	else if (xhtmlData && !(xhtmlData.indexOf('<span') >= 0)) {
		jQuery( "#patientinform-links").replaceWith(xhtmlData);
		fixColHeights(1);
		if (cbA.length) {
			updateCBItem(cbA, '<a href="#patientinform-links"><em>patient</em>INFORMation</a>', true);
		}
	}
	else {
		if (cbA.length) {
			updateCBItem(cbA, '<div id="cb-loaded-patientinform-none">No <em>patient</em>INFORMation available for this article</div>', false);
		}
	}
}
function addHWCiting(xhtmlData) {
	var cbA = jQuery( "#cb-loading-hw-cited");
	if (gIsFrameset) {
		if (cbA.length) {
			updateCBItem(cbA, '<div id="cb-loaded-hw-cited-none">Not available in this view</div>', false);
		}
	}
	else if (xhtmlData) {
		jQuery( "#panel-content-block").append(xhtmlData);
		var hwCitingLabel = getSiteOption('hwCitingLabel', 'View citing article information');

		fixColHeights(1);
		if (cbA.length) {
            if (!(getSiteOption("includeHWCitingTitle") == true)) {
                updateCBItem(cbA, '<a href="#cited-by">' + hwCitingLabel + '</a>', true);
            }
            else {
                updateCBItem(cbA, '<a href="#cited-by" title="HighWire Press-hosted articles citing this article">' + hwCitingLabel + '</a>', true);

            }
		}
	}
	else {
if (cbA.length) {
		updateCBItem(cbA, '<div id="cb-loaded-hw-cited-none">No citing articles</div>', false);
	    jQuery( "#cb-art-cited").addClass("no-citations");
	}
	}
}


function addISIRelated(xhtmlData) {
	var cbA = jQuery( "#cb-isi-similar-articles");
	if (xhtmlData && !(xhtmlData.indexOf('<span') >= 0)) {
		if (cbA.length) {
			updateCBItem(cbA, xhtmlData, true);
		}
	}
	else {
		if (cbA.length) {
			updateCBItem(cbA, '<div id="cb-loaded-isi-related-none">No Web of Science related articles</div>', false);
		}
	}
}
function addISICiting(xhtmlData) {
	var cbA = jQuery( "#cb-loading-isi-cited");
	if (xhtmlData && !(xhtmlData.indexOf('<span') >= 0)) {
		if (cbA.length) {
			updateCBItem(cbA, xhtmlData, true);
		}
	}
	else {
		if (cbA.length) {
			updateCBItem(cbA, '<div id="cb-loaded-isi-cited-none">No Web of Science citing articles</div>', false);
		}
	}
}
function addScopusCiting(xhtmlData) {
	var cbA = jQuery( "#cb-loading-scopus-cited");
	if (xhtmlData && (xhtmlData.indexOf('<a ') >= 0)) {
		if (cbA.length) {
			updateCBItem(cbA, xhtmlData, true);
		}
	}
	else {
		if (cbA.length) {
			updateCBItem(cbA, '<div id="cb-loaded-scopus-cited-none">No Scopus citing articles</div>', false);
		}
	}
}
function addEntrezLinks(xhtmlData) {
	var entrezDiv = jQuery( "#cb-entrez-links-placeholder");

	if (xhtmlData && (xhtmlData.indexOf('<a ') >= 0)) {
		if (entrezDiv.length) {
			updateCBItem(entrezDiv, xhtmlData, true);
			jQuery( entrezDiv).parent('li').addClass('has-data');
		}
	}
	else {
		if (entrezDiv.length) {
			updateCBItem(entrezDiv, '<div id="cb-entrez-links-none">No NCBI links</div>', false);
		}
	}
}

function updateSBLinks() {

	var fbLink = jQuery( "a.sb-facebook");
	if (fbLink.length) {
		fbLink.click(
			function(e) {
				window.open(this.href, 'sharer', 'toolbar=0,status=0,width=626,height=436');
				e.preventDefault();
			}
		);
	}

    updateGPlus();

}

function updateGPlus() {

    var  gPlus = jQuery( "ul.social-bookmark-links li.social-bookmarking-item-googleplus");
    /*
                See instructions here if you want to see how to change this:
                http://www.google.com.insb.bib.cnrs.fr/webmasters/+1/button/
            */
    /*
    	running a test here to check for MSIE < v8, in the interest of avoiding odd display when code returned by google fails.
    */
    var suppressGPlus = false;
    if (($.browser.msie)&&(parseInt($.browser.version, 10)<8)) {
    	suppressGPlus = true;
    }
    if (suppressGPlus) {
    	jQuery( '.social-bookmarking-item-googleplus').hide();
    } else {
		if (gPlus.length) {
	        var gPlusSize         = getSiteOption("googlePlusSize", 'small');
	        var gPlusDisplayCount = getSiteOption("googlePlusDisplayCount", 'false');
	        var gPlusURL          = jQuery( 'head meta[name=citation_public_url]').filter(':first').attr('content') || document.location;

	        //gPlus.prepend('<g:plusone size="' + gPlusSize + '" count="' + gPlusDisplayCount + '" href="' + gPlusURL +'" callback="gPlusCallback"></g:plusone>');
	        gPlus.prepend('<div class="g-plusone" data-size="' + gPlusSize + '" data-count="' + gPlusDisplayCount + '" data-href="' + gPlusURL +'" data-callback="gPlusCallback"></div>');

	        jQuery( 'a', gPlus).hide();  // remove anchor tag, we'll need it later for clicktracking

	        jQuery( 'body').append('<script type="text/javascript" src="https://apis.google.com/js/plusone.js"></script>');
		};
	}
}

function gPlusCallback() {
    var gPlusLoggerURL = jQuery( "ul.social-bookmark-links li.social-bookmarking-item-googleplus a").filter(':first').attr('href');

   // silently log the successful click
   if (gPlusLoggerURL.length) {
    $.get(gPlusLoggerURL);
   }

}

function addDockedNav() {
	var newDiv = '<div id="docked-nav"></div>';
	jQuery( "#col-2").append(newDiv);
	jQuery( "#col-2 #docked-nav").hide();
	var newDivJQuery = jQuery( "#docked-nav");
	newDivJQuery.append(jQuery( ".pane-highwire-article-views").clone());
	newDivJQuery.append(jQuery( ".pane-highwire-article-citation").clone());
	newDivJQuery.append(jQuery( ".pane-highwire-article-nav").clone());
	jQuery( "#col-2 #docked-nav").fadeIn(250);
}


function removeDockedNav() {
	var dockedNav = jQuery( "div#docked-nav");
	if(dockedNav.length) {
		dockedNav.fadeOut(250, function() { dockedNav.remove(); });
	}
}

function fixWackyReflinksMarkup() {
    // There's whitespace between the <li> tags, need to remove to avoid ugly spaces between words.


    jQuery( "div.ref-list ol.cit-auth-list,div.ref-list ol.cit-ed-list").each(
		function(i) {
			var original_html = jQuery( this).html();
			var whitespace_stripped_html;

            //Multiple spaces into one, remove trailing spaces after </li>
            whitespace_stripped_html = original_html.replace(/\s+/g,' ');
            whitespace_stripped_html = whitespace_stripped_html.replace(/<\/li>\s+/gi,'</li>');
            whitespace_stripped_html = whitespace_stripped_html.replace(/<\/span>\s+<\/li>/gi,'</span></li>');
			jQuery( this).html(whitespace_stripped_html);
		}
	);

}

function linkPDFExtImg() {
	var pdfExtImg = jQuery( "#panel-content-block div.extract-view img.pdf-extract-img");
	if (pdfExtImg.length) {
		pdfExtImg.before(
			'<p class="pdf-extract-click-text">' +
			getSiteOption('pdfExtractExpandText','Click image below to view at full size.') +
			'<\/p>'
		);
		pdfExtImg.wrap('<a class="pdf-extract-click" href="#">');
		var clickA = jQuery( "#panel-content-block div.extract-view a.pdf-extract-click");
		clickA.click(
			function(e) {
				var wasExpanded = jQuery( this).hasClass("expanded");
				if (wasExpanded) {
					jQuery( this).removeClass("expanded");
				}
				else {
					jQuery( this).addClass("expanded");
				}
				jQuery( "#content-option-box #content-toggle a").trigger("click");
				jQuery( "#panel-content-block div.extract-view p.pdf-extract-click-text").text(

					(wasExpanded ? getSiteOption('pdfExtractExpandText','Click image below to view at full size.') : getSiteOption('pdfExtractContractText','Click image below to return to normal size.'))
				);
				e.preventDefault();
			}
		);
	}
}
function AuthClickOnAffilButton(authorAffilList,authAffilMatch,expandStr,contractStr,e)
{
	var allViewMores = jQuery( authorAffilList);
	var authAffils = jQuery( authAffilMatch);
	if ((jQuery( authorAffilList).filter(':first')).text() == contractStr) {
		/* hide the affil list */
		allViewMores.empty().append(expandStr);
		authAffils.each(
			function(i) {
				modClass(authAffils.eq(i),'hideaffil','showaffil');
			}
		);
	}
	else {
		allViewMores.empty().append(contractStr);
		authAffils.each(
			function(i) {
				modClass(authAffils.eq(i),'showaffil','hideaffil');
			}
		);
	}
	fixColHeights(1);
	e.preventDefault();
}

function AuthClickOnNotesButton(authorNotes,authNotesMatch,expandStr,contractStr,e)
{
	var allViewMores = jQuery( authorNotes);
	var authNotes = jQuery( authNotesMatch);
	if ((jQuery( authorNotes).filter(':first')).text() == contractStr) {
		/* hide the affil list */
		allViewMores.empty().append(expandStr);
		authNotes.each(
			function(i) {
				modClass(authNotes.eq(i),'hidenotes','shownotes');
			}
		);
	}
	else {
		allViewMores.empty().append(contractStr);
		authNotes.each(
			function(i) {
				modClass(authNotes.eq(i),'shownotes','hidenotes');
			}
		);
	}
	fixColHeights(1);
	e.preventDefault();
}

function AuthClickOnContribAffilLink(authorAffilList,authAffilMatch,expandStr,contractStr,e)
{
	jQuery( authorAffilList).each(
		function() {
			if (jQuery( this).text() == expandStr) {
				jQuery( this).empty().append(contractStr);
				var authAffils = jQuery( authAffilMatch);
				if (authAffils.length) {
					authAffils.each(
						function(i) {
							modClass(authAffils.eq(i),'showaffil','hideaffil');
						}
					);
				}
				fixColHeights(1);
			}
		}
	);
}


function linkNewWin(jQEl) {
	var id = jQEl.attr("id");
	var targetName = "_blank";
	var config = null;
	if ((id != undefined) && id && !(id == '') &&
		!(gSiteOptions.openWindowDetails == undefined)) {
		var newLinkOptions = gSiteOptions.openWindowDetails[id];
		if ((newLinkOptions != undefined) && newLinkOptions) {
			var overrideTarget = newLinkOptions.target;
			var overrideConfig = newLinkOptions.config;
			if ((overrideTarget != undefined) && overrideTarget && (overrideTarget.length > 0)) {
				targetName = overrideTarget;
			}
			if ((overrideConfig != undefined) && overrideConfig && (overrideConfig.length > 0)) {
				config = overrideConfig;
			}
		}
	}
	var origTitle = jQEl.attr("title");
	var newTitle = '';
	if ((origTitle == undefined) || (!origTitle)) {
		origTitle = '';
	}
	else {
		newTitle = origTitle + ' ';
	}
	newTitle += '[opens in a new window]';
	jQEl.attr("target", targetName).attr("title", newTitle);
	if (config != null) {
		jQEl.click(function() {
			window.open(jQEl.attr("href"), targetName, config);
			return false;
		});
	}
}

function setupCollapsibles() {
	// sites can override this if they want to target different elements
	prepCollapsibles("div.content-box div.collapsible");
}

/* handle collapsing content-box areas, with optional user prefs */
function prepCollapsibles(parentJQueryEl, setPrefs) {
	if (typeof(setPrefs) == "undefined") {
		setPrefs = true;
	}
	var expandStr, contractStr;
	if (gSiteOptions == undefined) {
		expandStr = '+';
		contractStr = '-';
	}
	else {
		if (gSiteOptions.cbExpandString != undefined) {
			expandStr = gSiteOptions.cbExpandString;
		}
		else if (gSiteOptions.expandString != undefined) {
			expandStr = gSiteOptions.expandString;
		}
		else { expandStr = '+'; }
		if (gSiteOptions.cbContractString != undefined) {
			contractStr = gSiteOptions.cbContractString;
		}
		else if (gSiteOptions.contractString != undefined) {
			contractStr = gSiteOptions.contractString;
		}
		else { contractStr = '-'; }
	}
	jQuery( parentJQueryEl).each(
		function(i) {
			var $this = jQuery( this);
			var id = $this.attr("id");
			if (id != null) {
				var $labelElement, $label;
				for (var i = 0; i < gSiteOptions.collapsibleLabels.length; i++) {
					$labelElement = gSiteOptions.collapsibleLabels[i];
					$label = $this.find($labelElement);
					if ($label.length) {
						// get first one
						if ($label.length > 1) {
							$label = $label.eq(0);
						}
						break;
					}
				}

				// if there's already a toggle leave this as is
				var alreadyHasToggle = ($label.length && $label.find("a.collapse-toggle").length);
				if ($label.length) {
					if (!(alreadyHasToggle)) {
						var labelText = $label.html();
						$label.empty().append('<a href="#" class="collapse-toggle"><span class="view-more">' + contractStr + '</span> ' + labelText + '</a>');
					}
				}

				var isOpen = true;
				if (setPrefs) {
					var statePref = getPrefValue(id + "-state");
					if (statePref && ((statePref == 'open') || (statePref == 'closed'))) {
						if (statePref == 'closed') {
							isOpen = false;
						}
					}
					else {
						// no pref, or unknown value
						if ($this.hasClass('default-closed')) {
							isOpen = false;
						}
					}
				}
				var childListType;
				if ($this.children("ol").length) {
					childListType = "ol";
				}
				else {
					childListType = "ul";
				}
				if (!isOpen && !(alreadyHasToggle)) {
					toggleCollapse($this, childListType, ($labelElement + " a span.view-more"), expandStr, contractStr);
				}
				$this.find("a.collapse-toggle").click(
					function(e) {
						// need to find first parent that is in either state...
//						jQuery( this).parents(".collapsible, .collapsed").each(
						jQuery( this).parents().each(
							function(i) {
								if (jQuery( this).hasClass('collapsible') || jQuery( this).hasClass('collapsed')) {
									var isCollapsed = $this.hasClass('collapsed');
									if (setPrefs) {
										var prefName = id + "-state";
										if (isCollapsed) {
											setPref(prefName, "open");
										}
										else {
											setPref(prefName, "closed");
										}
									}
									toggleCollapse(jQuery( this), childListType, ($labelElement + " a span.view-more"), expandStr, contractStr);
									fixColHeights(-1);
									return false;
								}
							}
						);
						e.preventDefault();
					}
				);
			}
		}
	);
}

function handleMostBoxes() {
	var hadMostBoxes = false;
	jQuery( "div.most-links-box").each(
		function(i) {
			var $this = jQuery( this);
			modClass($this, '', 'js-marker');
			hadMostBoxes = true;
			// 1) get the <ul> child
			var $mostlist = $this.children("ul");
			// 2) build a new <ul> for the header, from <h4> child:
			var hasSelected = false;
			var mostHdrList = null;
			if ($mostlist.length) {
				var hasItems = false;
				var hdrList = '<ul class="most-headings">';
				var $mostli = $mostlist.children("li");
				if ($mostli.length) {
					hasItems = true;
					$mostli.each(
						function(i) {
							var $hdr = jQuery( this).children("h4").html();
							if (jQuery( this).hasClass("most-cur-sel")) {
								hasSelected = true;
								hdrList += '<li class="' + jQuery( this).attr("class") + '"><a href="#">' + $hdr + '</a></li>';
							}
							else {
								hdrList += '<li><a href="#">' + $hdr + '</a></li>';
							}
						}
					);
				}
				hdrList += '</ul>';
				if (hasItems) {
					var header = $this.find("div.most-header h3");
					if (header.length) {
						header.after(hdrList);
					}
				}
				mostHdrList = $this.find("div.most-header ul");
			}
			// 3) if none of the items are selected, make the first one be
			if (((mostHdrList != null) && mostHdrList.length) && (!hasSelected)) {
				modClass(mostHdrList.find("li:first"), "most-cur-sel", "");
				modClass($this.children("ul").children("li:first"), "most-cur-sel", "");
			}
			// 4) bind a handler to each heading li a to toggle class="most-cur-sel" for
			//    that item and its corresponding item in the list below
			if (((mostHdrList != null) && mostHdrList.length)) {
				mostHdrList.find("li a").click(
					function(e) {
						var curA = jQuery( this);
						toggleMostSelection(e, $this, curA.text());
					}
				);
			}
		}
	);
	if (hadMostBoxes) {
		fixColHeights(-1);
	}
}

function toggleMostSelection(e, mostItemDiv, groupName) {
	if (mostItemDiv.length) {
		var hdrItems = mostItemDiv.find("div.most-header ul li");
		if (hdrItems.length) {
			hdrItems.each(
				function() {
					var hdrName = jQuery( this).text();
					var isChoice = (hdrName == groupName);
					if (isChoice) {
						modClass(jQuery( this), 'most-cur-sel', '');
					}
					else {
						modClass(jQuery( this), '', 'most-cur-sel');
					}
				}
			);
		}
		var listItems = mostItemDiv.children("ul").find("li");
		if (listItems.length) {
			listItems.each(
				function() {
					var itemName = jQuery( this).children("h4").text();
					var isChoice = (itemName == groupName);
					if (isChoice) {
						modClass(jQuery( this), 'most-cur-sel', '');
					}
					else {
						modClass(jQuery( this), '', 'most-cur-sel');
					}
				}
			);
		}
	}
	fixColHeights(-1);
	e.preventDefault();
}

function toggleCollapse(collapsibleEl, hideTagExpr, modTextTagExpr, expandStr, contractStr) {
	if (typeof(expandStr) == "undefined") {
		expandStr = getSiteOption('expandString','+');
	}
	if (typeof(contractStr) == "undefined") {
		contractStr = getSiteOption('contractString','-');
	}

	var $this = collapsibleEl;
	if ($this.hasClass('collapsed')) {
		modClass($this, 'collapsible', 'collapsed');
	}
	else if ($this.hasClass('collapsible')) {
		modClass($this, 'collapsed', 'collapsible');
	}
	var $toToggle = $this.find(hideTagExpr);
	if ($toToggle.length) {
		$toToggle.toggle();
	}
	var $rewrite = $this.find(modTextTagExpr);
	if ($rewrite.length) {
		$rewrite.each(
			function(i, e) {
				var txt = jQuery( e).text();
				if (txt == expandStr) {
					jQuery( e).empty().append(contractStr);
				}
				else if (txt == contractStr) {
					jQuery( e).empty().append(expandStr);
				}
			}
		);
	}
}

function setupDockBlock(col, dockId, dockClass, contentRules) {
	if ((col == 2) || (col == 3)) {
		var navOrigHeight = ((col == 2) ? gColOrigHeights.col2 : gColOrigHeights.col3);
		var gColNum = col;
		var gNav = ((col == 2) ? gCol2 : gCol3);
		var gNavHeight = ((gNav.length) ? navOrigHeight : 0);
		var gNavTopOffset =  (gNavHeight > 0 ? gNav.offset().top : 0);
		var gScrollPos = 0;
		var gNavDocked = false;
		/* this won't work in IE 6, so skip it */
		/* article navigation docking block */
		jQuery( window).scroll(function () {
			if (gNavHeight > 0 && !isIE6) {
				var offsets = getPageOffset();
				var origHeight = ((gColNum == 2) ? gColOrigHeights.col2 : gColOrigHeights.col3);
				var navBottom = (gNavTopOffset + origHeight) - offsets.y;
				if (gScrollPos != offsets.y) {
					gScrollPos = offsets.y;
					if ((navBottom <= 0) && !gNavDocked) {
						/* dock the nav */
						gNavDocked = true;
						addDockBlock(col, dockId, dockClass, contentRules);
						addDockBlockCallback(dockId);
					}
					else if (gNavDocked && (navBottom > 0)) {
						/* undock the nav */
						gNavDocked = false;
						removeDockBlock(col, dockId);
						removeDockBlockCallback(dockId);
					}
				}
			}
		});
	}
}
function addDockBlock(col, dockId, dockClass, contentRules) {
	// get direct children li elements only
	var colId = ("col-" + col);
	var newDiv = '<div' + ((dockId == '') ? '': (' id="' + dockId + '"')) + ((dockClass == '') ? '': (' class="' + dockClass + '"')) + '></div>';
	jQuery( "#" + colId).append(newDiv);
	jQuery( "#" + colId + " #" + dockId).hide();
	var newDivJQuery = jQuery( "#" + dockId);
	if (contentRules.length > 1) {
		for (var i = 0; i < (contentRules.length - 1); i = i + 2) {
			var addTo = processDockingRule(contentRules[i]);
			var add = processDockingRule(contentRules[i + 1]);
			var addTarget, addVal;
			if (addTo.ruleType == 'empty') {
				addTarget = newDivJQuery;
			}
			else if (addTo.ruleType == 'jQuery') {
				addTarget = jQuery( addTo.rule);
			}
			if (add.ruleType = 'jQuery') {
				var addJQEl = jQuery( add.rule);
				if (addJQEl.length) {
					addVal = addJQEl.clone();
				}
				else {
					addVal = '';
				}
			}
			else if (add.ruleType = 'html') {
				addVal = add.rule;
			}
			if (addTarget.length && !(addVal == '')) {
				addTarget.append(addVal);
			}
			else {
			}
		}
	}
	var newDockJQrule = "#" + colId + " #" + dockId;
	var newDock = jQuery( newDockJQrule);
	if (newDock.length) {
		prepCollapsibles(newDockJQrule + " div.collapsible, " + newDockJQrule + " div.collapsed", false);
		newDock.fadeIn(250);
	}
}

function processDockingRule(rule) {
	var retObj = new Object();
	retObj.empty = (((rule != undefined) && (!(rule == ''))) ? false : true);
	if (rule == '') {
		retObj.ruleType = 'empty';
	}
	else if (/^\$\(.+\)$/.test(rule)) {
		retObj.ruleType = 'jQuery';
		retObj.rule = rule.replace(/^\$\((.+)\)$/,"$1");
	}
	else {
		retObj.ruleType = 'html';
		retObj.rule = rule;
	}
	return retObj;
}

function removeDockBlock(col, dockId) {
	var dockedNav = jQuery( "div#" + dockId);
	if(dockedNav.length) {
		dockedNav.fadeOut(250, function() { dockedNav.remove(); });
	}
}
// custom implementations can override this function to do something when a dockBlock is added
// they should return false if they want to stop any additional default handling
function customAddDockBlockCallback(dockId) {
	return true;
}
// this function should not be overridden
function addDockBlockCallback(dockId) {
	if (customAddDockBlockCallback(dockId) != false) {
		// no default handling, currently
	}
}
// custom implementations can override this function to do something when a dockBlock is removed
// they should return false if they want to stop any additional default handling
function customRemoveDockBlockCallback(dockId) {
	return true;
}
// this function should not be overridden
function removeDockBlockCallback(dockId) {
	if (customRemoveDockBlockCallback(dockId) != false) {
		// no default handling, currently
	}
}

function updateFormInput(labelMatchString, inputMatchString, defaultColorString, textColorString) {
	if ((defaultColorString == null) || (defaultColorString == '')) {
		defaultColorString = "#A0A0A0";
	}
	if ((textColorString == null) || (textColorString == '')) {
		textColorString = "black";
	}
	var label = jQuery( labelMatchString);
	var input = jQuery( inputMatchString);
	if (input.length) {
		if ((label.length) && (input.val() == '')) {
			input.val(label.text()).css("color",defaultColorString);
		}
		input.focus(
			function(e) {
				if ((label.length) &&  (jQuery( this).val() == label.text())) {
					jQuery( this).val('').css("color",textColorString);
				}
			}
		);
		input.blur(
			function(e) {
				if ((label.length) && (jQuery( this).val() == '')) {
					jQuery( this).val(label.text()).css("color",defaultColorString);
				}
			}
		);
	}
	var parentForm = label.parents("form").eq(0);
	if (parentForm.length) {
		parentForm.submit(
			function() {
				if ((label.length) && (input.length) && (input.val() == label.text())) {
					input.val('');
				}
				return true;
			}
		);
	}
}

function checkUnloadedImgs(lookupRule) {
	if (lookupRule == null) {
		lookupRule = "img";
	}
	var imageNum = -1;
	var images = jQuery( lookupRule);
	if (images.length) {
		imageNum = 0;
		images.each(
			function() {
				//if (jQuery( this).height() == 0) {
				if (this.offsetHeight == 0) {
					imageNum++;
				}
			}
		);
	}
	return imageNum;
}

function fixHeightForImages(iter, numImagesToLoad, lookupRule) {
	if (lookupRule == null) {
		lookupRule = "img";
	}
	if (numImagesToLoad > 0) {
		var lastNumImagesToLoad = numImagesToLoad;
		var newNumImagesToLoad = checkUnloadedImgs(lookupRule);
		if (newNumImagesToLoad < lastNumImagesToLoad) {
			fixColHeights(1);
			if ((newNumImagesToLoad > 0) && (iter < 10)) {
				setTimeout("fixHeightForImages(" + (iter+1) + "," + newNumImagesToLoad + ",'" + lookupRule + "')", 1000);
			}
		}
	}
}
function getColInfo(colJQEl) {
	var colInfo = { valid: false, height: 0, bottom: 0, lastElBottom: 0, extra: 0, fix: false };
	if (colJQEl.length) {
		var c = colJQEl.children();
		if (c.length) {
			colInfo.height = getObjHeight(colJQEl);
			colInfo.bottom = colJQEl.offset().top + colInfo.height;
			var last = c.eq(c.length - 1);
			colInfo.lastElBottom = last.offset().top + getObjHeight(last);
			colInfo.extra = colInfo.bottom - colInfo.lastElBottom;
			colInfo.valid = true;
		}
	}
	return colInfo;
}
function checkColHeights() {
	var colToFix = 0;
	if ((gIsMultiCol || gColTempResize) && ((gColOrigHeights.main + gColOrigHeights.col2 + gColOrigHeights.col3) > 0)) {
		var colMainInfo = getColInfo(gColMain);
		var col2Info = getColInfo(gCol2);
		var col3Info = getColInfo(gCol3);
		if (colMainInfo.valid && (gColLows.main == 0)) { gColLows.main = colMainInfo.lastElBottom; }
		if (col2Info.valid && (gColLows.col2 == 0)) { gColLows.col2 = col2Info.lastElBottom; }
		if (col3Info.valid && (gColLows.col3 == 0)) { gColLows.col3 = col3Info.lastElBottom; }
		// check if low point in any column has changed
		if ((colMainInfo.valid) && (colMainInfo.lastElBottom != gColLows.main)) {
			colMainInfo.fix = true;
		}
		if ((col2Info.valid) && (col2Info.lastElBottom != gColLows.col2)) {
			col2Info.fix = true;
		}
		if ((col3Info.valid) && (col3Info.lastElBottom != gColLows.col3)) {
			col3Info.fix = true;
		}
		// check for columns not having same bottom points
		if ((colMainInfo.valid) &&
			  (((col2Info.valid) && (colMainInfo.height < col2Info.height)) ||
			   ((col3Info.valid) && (colMainInfo.height < col3Info.height)))) {
			colMainInfo.fix = true;
		}
		if ((col2Info.valid) &&
			  (((colMainInfo.valid) && (col2Info.height < colMainInfo.height)) ||
			   ((col3Info.valid) && (col2Info.height < col3Info.height)))) {
			col2Info.fix = true;
		}
		if ((col3Info.valid) &&
			  (((colMainInfo.valid) && (col3Info.height < colMainInfo.height)) ||
			   ((col2Info.valid) && (col3Info.height < col2Info.height)))) {
			col3Info.fix = true;
		}

		if (colToFix > -1) {
			if (colMainInfo.valid && colMainInfo.fix) {
				if ((col2Info.valid && col2Info.fix) || (col3Info.valid && col3Info.fix)) {
					colToFix = -1;
				}
				else { colToFix = 1; }
			}
			else if (col2Info.valid && col2Info.fix) {
				if (col3Info.valid && col3Info.fix) {
					colToFix = -1;
				}
				else { colToFix = 2; }
			}
			else if (col3Info.valid && col3Info.fix) {
				colToFix = 3;
			}
		}

		if (colToFix != 0) {
			fixColHeights(colToFix);
		}
		if (colMainInfo.valid) { gColLows.main = colMainInfo.lastElBottom; }
		if (col2Info.valid) { gColLows.col2 = col2Info.lastElBottom; }
		if (col3Info.valid) { gColLows.col3 = col3Info.lastElBottom; }
	}
	setTimeout("checkColHeights()", 5000);
}
function fixColHeights(colChanged) {
	if ((gIsMultiCol || gColTempResize) && getSiteOption('setColSizes', true)) {
		var setAll = false;
		if (colChanged == 1) {
			if (gColMain.length) {
				setObjHeight(gColMain,'auto');
				gColOrigHeights.main = getObjHeight(gColMain);
			}
		}
		else if (colChanged == 2) {
			if (gCol2.length) {
				setObjHeight(gCol2,'auto');
				gColOrigHeights.col2 = getObjHeight(gCol2);
			}
		}
		else if (colChanged == 3) {
			if (gCol3.length) {
				setObjHeight(gCol3,'auto');
				gColOrigHeights.col3 = getObjHeight(gCol3);
			}
		}
		else {
			if (colChanged == -1) {
				// force resize all 3
				if (gColMain.length) { setObjHeight(gColMain,'auto'); }
				if (gCol2.length) { setObjHeight(gCol2,'auto'); }
				if (gCol3.length) { setObjHeight(gCol3,'auto'); }
			}
			gColOrigHeights.main = ((gColMain.length) ? (getObjHeight(gColMain)) : 0);
			gColOrigHeights.col2 = ((gCol2.length) ? (getObjHeight(gCol2)) : 0);
			gColOrigHeights.col3 = ((gCol3.length) ? (getObjHeight(gCol3)) : 0);
			setAll = true;
		}
		if ((gColOrigHeights.main + gColOrigHeights.col2 + gColOrigHeights.col3) > 0) {
			var colInfo;
			var maxH = Math.max(Math.max(gColOrigHeights.main, gColOrigHeights.col2), gColOrigHeights.col3);
			if ((gColMain.length) && ((gColOrigHeights.main < maxH) || (getObjHeight(gColMain) > maxH) || setAll)) {
				setObjHeight(gColMain, maxH);
				colInfo = getColInfo(gColMain);
				if (colInfo.valid) { gColLows.main = colInfo.lastElBottom; }
			}
			if ((gCol2.length) && ((gColOrigHeights.col2 < maxH) || (getObjHeight(gCol2) > maxH) || setAll)) {
				setObjHeight(gCol2, maxH);
				colInfo = getColInfo(gCol2);
				if (colInfo.valid) { gColLows.col2 = colInfo.lastElBottom; }
			}
			if ((gCol3.length) && ((gColOrigHeights.col3 < maxH) || (getObjHeight(gCol3) > maxH) || setAll)) {
				setObjHeight(gCol3, maxH);
				colInfo = getColInfo(gCol3);
				if (colInfo.valid) { gColLows.col3 = colInfo.lastElBottom; }
			}
		}
	}
}

function getObjHeight(obj) {
	if (obj instanceof jQuery) {
		return obj.get(0).offsetHeight;
	}
	else if (obj instanceof Array) {
		return obj[0].offsetHeight;
	}
	else {
		return obj.offsetHeight;
	}
}
function getObjWidth(obj) {
	if (obj instanceof jQuery) {
		return obj.get(0).offsetWidth;
	}
	else if (obj instanceof Array) {
		return obj[0].offsetWidth;
	}
	else {
		return obj.offsetWidth;
	}
}
function setObjHeight(obj, h) {
	if (typeof(h) == 'number') {
		h = "" + h + "px";
	}
	if (obj instanceof jQuery) {
		obj.get(0).style.height = h;
	}
	else if (obj instanceof Array) {
		obj[0].style.height = h;
	}
	else {
		obj.style.height = h;
	}
}
function setObjWidth(obj, w) {
	if (typeof(w) == 'number') {
		w = "" + w + "px";
	}
	if (obj instanceof jQuery) {
		obj.get(0).style.width = w;
	}
	else if (obj instanceof Array) {
		obj[0].style.width = w;
	}
	else {
		obj.style.width = w;
	}
}

/* num pixels top of page has been scrolled offscreen */
function getPageOffset() {
	var offset = new Object();
  if( typeof( window.pageXOffset ) == 'number' ) {
    offset.x = window.pageXOffset;
    offset.y = window.pageYOffset;
  } else if( document.body && ( document.body.scrollTop ) ) {
    offset.x = document.body.scrollLeft;
    offset.y = document.body.scrollTop;
  } else if( document.documentElement && ( document.documentElement.scrollTop ) ) {
    offset.x = document.documentElement.scrollLeft;
    offset.y = document.documentElement.scrollTop;
  }
  if (typeof(offset.x)=="undefined") { offset.x = 0; }
  if (typeof(offset.y)=="undefined") { offset.y = 0; }
  return offset;
}
function getViewportDim() {
	var dim = new Object();
	// non-IE
	if (typeof window.innerWidth != 'undefined') {
		dim.x = window.innerWidth;
		dim.y = window.innerHeight;
	}
	// IE6 in standards compliant mode
	else if (typeof document.documentElement != 'undefined'
		&& typeof document.documentElement.clientWidth !=
		'undefined' && document.documentElement.clientWidth != 0) {
		dim.x = document.documentElement.clientWidth;
		dim.y = document.documentElement.clientHeight;
	}
	// older versions of IE
	else {
		dim.x = document.getElementsByTagName('body')[0].clientWidth;
		dim.y = document.getElementsByTagName('body')[0].clientHeight;
	}
	return dim;
}

function newWindowTargets() {
	var newWins = jQuery( "a.in-nw");
	if (newWins.length) {
		newWins.each(
			function() {
				var $this = jQuery( this);
				$this.attr("target", "_blank");
				modClass($this, 'in-nw-vis', 'in-nw');
			}
		);
	}
}

function getSiteOption(optionName, defaultVal) {
	if (typeof(defaultVal) == "undefined") {
		defaultVal = gSiteOptionUnknown;
	}
	if (gSiteOptions == undefined) {
		return defaultVal;
	}
	else if (gSiteOptions[optionName] == undefined) {
		return defaultVal;
	}
	else {
		return gSiteOptions[optionName];
	}
}

function addPartHeaders(req) {
	req.setRequestHeader('Accept', 'application/xhtml+xml');
	req.setRequestHeader('Range', 'part=variant-contents');
}
function addCommonHeaders(req) {
	if (typeof(callbackToken) != 'undefined') {
		req.setRequestHeader('X-Token', callbackToken);
	}
}
function allowsCookies() {
	if (hasTestCookie()) {
		return true;
	}
	else {
		setTestCookie();
		return hasTestCookie();
	}
}
function setTestCookie() {
	setCookie('cks','allowed',7,'/',null,null);
}
function hasTestCookie() {
	var ckVal = getCookie('cks');
	if ((ckVal != null) && (ckVal == 'allowed')) {
		return true;
	}
	else {
		return false;
	}
}

var UIPrefsCk={
	name: "UIPrefs",
	expDays: 3652,
	path: '/'
};

function prefDefined(name) {
	var val = getPrefValue(name);
	if (val == null) {
		return false;
	}
	else {
		return true;
	}
}
function getPrefValue(name) {
	var prefArray = getPrefArray();
	if (prefArray) {
		name = convertPrefString(name);
		var match = name + ":";
		for (var i = 0; i < prefArray.length; i++) {
			if (prefArray[i].indexOf(match) == 0) {
				return prefArray[i].substring(match.length);
			}
			else if (prefArray[i] == name) {
				// name only, no value, return boolean true
				return true;
			}
		}
	}
	return null;
}
function convertPrefString(str) {
	if (str) {
		return str.replace(/[,:=]/g,"");
	}
	return str;
}
function removePref(name) {
	var prefArray = getPrefArray();
	name = convertPrefString(name);
	if (prefDefined(name)) {
		var newArray = new Array(prefArray.length - 1);
		var match = name + ":";
		var newArrIndex = 0;
		for (var i = 0; i < prefArray.length; i++) {
			if (!((prefArray[i].indexOf(match) == 0) || (prefArray[i] == name))) {
				newArray[newArrIndex++] = prefArray[i];
			}
		}
		setPrefCookie(newArray);
		if (prefArray.length == 1) {
			return null;
		}
		else {
			return newArray;
		}
	}
	else {
		return prefArray;
	}
}

function setPref(name, value) {
	if (typeof(value) == "undefined") {
		value = null;
	}
	if (name) {
		name = convertPrefString(name);
		var prefArray = removePref(name);
		if (prefArray == null) {
			prefArray = new Array;
		}
		var newArray = new Array(prefArray.length + 1);
		var newPref = name + (((value != null) && (value != '')) ? (':' + value) : '');
		newArray[0] = newPref;
		for (var i = 0; i < prefArray.length; i++) {
			newArray[i+1] = prefArray[i];
		}
		setPrefCookie(newArray);
		return newArray;
	}
	else {
		return getPrefArray();
	}
}
function setPrefCookie(prefs) {
	if (prefs && (prefs.length > 0)) {
		var prefString;
		if (typeof(prefs) == 'string') {
			prefString = prefs;
		}
		else {
			if (prefs.length == 1) {
				prefString = prefs[0];
			}
			else {
				prefString = prefs.join(',');
			}
		}
		setCookie(UIPrefsCk.name, prefString, UIPrefsCk.expDays, UIPrefsCk.path, false, false);
	}
	else {
		deleteCookie(UIPrefsCk.name, UIPrefsCk.path, false);
	}
}
// individual pref items in the cookie are separated by ',' and
// nv pref pairs are joined by ':', i.e. n1:v1,n2:v2,n3:v3,etc.
// values are optional. (i.e. "n1,n2:v2,n3" is ok)
function getPrefArray() {
	var cookieVal = getCookie(UIPrefsCk.name);
	if (cookieVal != null) {
		return cookieVal.split(',');
	}
	return null;
}

// http://www.dustindiaz.com/top-ten-javascript/
function getCookie( name ) {
	var start = document.cookie.indexOf( name + "=" );
	var len = start + name.length + 1;
	if (((!start ) && (name != document.cookie.substring( 0, name.length ))) || (start == -1)) {
		return null;
	}
	var end = document.cookie.indexOf( ';', len );
	if ( end == -1 ) end = document.cookie.length;
	return unescape( document.cookie.substring( len, end ) );
}
function setCookie( name, value, expires, path, domain, secure ) {
	var today = new Date();
	today.setTime( today.getTime() );
	if ( expires ) {
		expires = expires * 1000 * 60 * 60 * 24;
	}
	var expires_date = new Date( today.getTime() + (expires) );
	document.cookie = name+'='+escape( value ) +
		( ( expires ) ? ';expires='+expires_date.toGMTString() : '' ) +
		( ( path ) ? ';path=' + path : '' ) +
		( ( domain ) ? ';domain=' + domain : '' ) +
		( ( secure ) ? ';secure' : '' );
}

function deleteCookie( name, path, domain ) {
	if ( getCookie( name ) ) document.cookie = name + '=' +
			( ( path ) ? ';path=' + path : '') +
			( ( domain ) ? ';domain=' + domain : '' ) +
			';expires=Thu, 01-Jan-1970 00:00:01 GMT';
}
function modClass(jQueryEl, addClass, removeClass) {
	if (jQueryEl.length) {
		var exClass = jQueryEl.attr("class");
		var classArr = new Array;
		if ((exClass != undefined) && (exClass != '')) {
			classArr = exClass.split(' ');
		}
		var newClass = '';
		var addClassExisted = false;
		for (var i = 0; i < classArr.length; i++) {
			var classPart = classArr[i];
			if (classPart != '') {
				if (classPart != removeClass) {
					if (classPart == addClass) {
						addClassExisted = true;
					}
					newClass = addToDelimString(newClass, classPart, ' ');
				}
			}
		}
		if (!addClassExisted) {
			newClass = addToDelimString(newClass, addClass, ' ');
		}
		jQueryEl.attr("class",newClass);
	}
}
function addToDelimString(origString, newPart, delim) {
	var newString;
	if ((origString == undefined) || (origString == '')) {
		newString = newPart;
	}
	else {
		newString = origString + delim + newPart;
	}
	return newString;
}
function debugOut(msg) {
	if(window.console) {
		window.console.log(msg);
	}
	/*
	else {
		alert(msg);
	}
	*/
}

    function is_defined( variable)
    {
        return (typeof(window[variable]) == "undefined")?  false: true;
    }

