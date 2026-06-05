function plot1 (AR_init, AR_noz, fpr, P_inf, P1, P2, P3, P4, P5,P6, T_inf,T1, T2, T3, T4, T5, T6, P_inf_tot, P1_0, P2_0, P3_0, P4_0, P5_0, P6_0, V_inf, v1, v2, v3, v4, v5, v6, v_channel_air, M_inf, M1, M2, M3, M4, M5, M6)


    x = [0, 1, 2, 3, 4, 5, 6];
    y1 = [P_inf, P1, P2, P3, P4, P5, P6];
    y2 = [P_inf_tot, P1_0, P2_0, P3_0, P4_0, P5_0, P6_0];
    y3 = [T_inf, T1, T2, T3, T4, T5, T6];
    x4 = [0, 1, 2, 3, 3.5, 4, 5, 6];
    y4 = [V_inf, v1, v2, v3, v_channel_air, v4, v5, v6];
    y5 = [M_inf, M1, M2, M3, M4, M5, M6];

    persistent ax1 ax2 ax3 ax4 ax5 fig created_count

    if isempty(created_count)
        
        fig = figure('Name','TMS Stations','NumberTitle','off');
        ax1 = subplot(5,1,1); hold(ax1,'on'); grid(ax1,'on');
        ax2 = subplot(5,1,2); hold(ax2,'on'); grid(ax2,'on');
        ax3 = subplot(5,1,3); hold(ax3,'on'); grid(ax3,'on');
        ax4 = subplot(5,1,4); hold(ax4,'on'); grid(ax4,'on');
        ax5 = subplot(5,1,5); hold(ax5,'on'); grid(ax5,'on');

        
        axs = [ax1 ax2 ax3 ax4 ax5];
        for a = axs
            a.FontSize = 12;
            a.XTick = x;
            xlabel(a,'Station','FontSize',10);
        end
        ylabel(ax1,'Static pressure','FontSize',10);
        ylabel(ax2,'Total pressure','FontSize',10);
        ylabel(ax3,'Static temperature','FontSize',10);
        ylabel(ax4,'Mach number','FontSize',10);
        ylabel(ax5,'Air velocity','FontSize',10);

        created_count = 0;
    end

    
    labelStr = sprintf('Diff. AR=%.2g, Noz. AR=%.2g', AR_init, AR_noz);

    
    plot(ax1, x, y1, '-o', 'LineWidth', 1.5, 'DisplayName', labelStr);
    plot(ax2, x, y2, '-o', 'LineWidth', 1.5, 'DisplayName', labelStr);
    plot(ax3, x, y3, '-o', 'LineWidth', 1.5, 'DisplayName', labelStr);
    plot(ax4, x, y5, '-o', 'LineWidth', 1.5, 'DisplayName', labelStr);
    plot(ax5, x4, y4, '-o', 'LineWidth', 1.5, 'DisplayName', labelStr);  

    created_count = created_count + 1;

   
    if created_count == 1
        legend(ax1,'Location','best');
        legend(ax2,'Location','best');
        legend(ax3,'Location','best');
        legend(ax4,'Location','best');
        legend(ax5,'Location','best');
    end

end