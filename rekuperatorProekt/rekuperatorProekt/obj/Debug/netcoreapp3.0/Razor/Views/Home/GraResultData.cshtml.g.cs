#pragma checksum "H:\proektRekuperator-master\rekuperatorProekt\rekuperatorProekt\Views\Home\GraResultData.cshtml" "{ff1816ec-aa5e-4d10-87f7-6f4963833460}" "6524566dc1fcf6af3d7c8540f9c74f8b4aca850f"
// <auto-generated/>
#pragma warning disable 1591
[assembly: global::Microsoft.AspNetCore.Razor.Hosting.RazorCompiledItemAttribute(typeof(AspNetCore.Views_Home_GraResultData), @"mvc.1.0.view", @"/Views/Home/GraResultData.cshtml")]
namespace AspNetCore
{
    #line hidden
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using System.Threading.Tasks;
    using Microsoft.AspNetCore.Mvc;
    using Microsoft.AspNetCore.Mvc.Rendering;
    using Microsoft.AspNetCore.Mvc.ViewFeatures;
#nullable restore
#line 1 "H:\proektRekuperator-master\rekuperatorProekt\rekuperatorProekt\Views\_ViewImports.cshtml"
using rekuperatorProekt;

#line default
#line hidden
#nullable disable
#nullable restore
#line 2 "H:\proektRekuperator-master\rekuperatorProekt\rekuperatorProekt\Views\_ViewImports.cshtml"
using rekuperatorProekt.Models;

#line default
#line hidden
#nullable disable
    [global::Microsoft.AspNetCore.Razor.Hosting.RazorSourceChecksumAttribute(@"SHA1", @"6524566dc1fcf6af3d7c8540f9c74f8b4aca850f", @"/Views/Home/GraResultData.cshtml")]
    [global::Microsoft.AspNetCore.Razor.Hosting.RazorSourceChecksumAttribute(@"SHA1", @"31d01f78bc8da9547bbe987991a0c300ca29fb43", @"/Views/_ViewImports.cshtml")]
    public class Views_Home_GraResultData : global::Microsoft.AspNetCore.Mvc.Razor.RazorPage<dynamic>
    {
        #pragma warning disable 1998
        public async override global::System.Threading.Tasks.Task ExecuteAsync()
        {
            WriteLiteral("\n");
#nullable restore
#line 2 "H:\proektRekuperator-master\rekuperatorProekt\rekuperatorProekt\Views\Home\GraResultData.cshtml"
  
    ViewData["Title"] = "Отчет";

#line default
#line hidden
#nullable disable
            WriteLiteral(@"
<div class=""text-center"">
    <h1 class=""display-4""><b>График</b></h1>

</div>

<canvas id=""myChart"" width=""400"" height=""150""></canvas>

<script>
    var ctx = document.getElementById('myChart');

    // Global Options:
    Chart.defaults.global.defaultFontColor = 'black';
    Chart.defaults.global.defaultFontSize = 16;

    var data = {
        labels: ");
#nullable restore
#line 21 "H:\proektRekuperator-master\rekuperatorProekt\rekuperatorProekt\Views\Home\GraResultData.cshtml"
           Write(Json.Serialize(ViewBag.x));

#line default
#line hidden
#nullable disable
            WriteLiteral(@",
        datasets: [{
            label: ""Stock A"",
            fill: false,
            lineTension: 0.1,
            backgroundColor: ""rgba(225,0,0,0.4)"",
            borderColor: ""red"", // The main line color
            borderCapStyle: 'square',
            borderDash: [], // try [5, 15] for instance
            borderDashOffset: 0.0,
            borderJoinStyle: 'miter',
            pointBorderColor: ""black"",
            pointBackgroundColor: ""white"",
            pointBorderWidth: 1,
            pointHoverRadius: 8,
            pointHoverBackgroundColor: ""yellow"",
            pointHoverBorderColor: ""brown"",
            pointHoverBorderWidth: 2,
            pointRadius: 4,
            pointHitRadius: 10,
            // notice the gap in the data and the spanGaps: true
            data:");
#nullable restore
#line 42 "H:\proektRekuperator-master\rekuperatorProekt\rekuperatorProekt\Views\Home\GraResultData.cshtml"
            Write(Json.Serialize(ViewBag.y));

#line default
#line hidden
#nullable disable
            WriteLiteral(@",
            spanGaps: true,
        }

        ]
    };

    // Notice the scaleLabel at the same level as Ticks
    var options = {
        scales: {
            yAxes: [{
                ticks: {
                    beginAtZero: true
                },
                scaleLabel: {
                    display: true,
                    labelString:    'Количество продуктов горения перед рекуператором',
                    fontSize: 20
                }
            }],

            xAxes: [{
                ticks: {
                    beginAtZero: true
                },
                scaleLabel: {
                    display: true,
                    labelString: 'Расход топлива на печь',
                    fontSize: 20
                }
            }]
        }
    };

    // Chart declaration:
    var myBarChart = new Chart(ctx, {
        type: 'line',
        data: data,
        options: options
    });
</script>");
        }
        #pragma warning restore 1998
        [global::Microsoft.AspNetCore.Mvc.Razor.Internal.RazorInjectAttribute]
        public global::Microsoft.AspNetCore.Mvc.ViewFeatures.IModelExpressionProvider ModelExpressionProvider { get; private set; }
        [global::Microsoft.AspNetCore.Mvc.Razor.Internal.RazorInjectAttribute]
        public global::Microsoft.AspNetCore.Mvc.IUrlHelper Url { get; private set; }
        [global::Microsoft.AspNetCore.Mvc.Razor.Internal.RazorInjectAttribute]
        public global::Microsoft.AspNetCore.Mvc.IViewComponentHelper Component { get; private set; }
        [global::Microsoft.AspNetCore.Mvc.Razor.Internal.RazorInjectAttribute]
        public global::Microsoft.AspNetCore.Mvc.Rendering.IJsonHelper Json { get; private set; }
        [global::Microsoft.AspNetCore.Mvc.Razor.Internal.RazorInjectAttribute]
        public global::Microsoft.AspNetCore.Mvc.Rendering.IHtmlHelper<dynamic> Html { get; private set; }
    }
}
#pragma warning restore 1591