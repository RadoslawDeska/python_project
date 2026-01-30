"""State Manager - Handles application state and flags.

Manages all state variables for measurement, fitting, and UI operations.
"""


class StateManager:
    """Manages application state variables and flags."""
    
    def __init__(self):
        """Initialize state manager with all application states."""
        # Measurement states
        self.clearing = False
        self.data_acquisition_complete = False
        self.experiment_stopped = False
        self.running = False
        
        # Initialization states
        self.initialized = False
        self.initializing = False

        # Data loading states
        self._file_loading = False
        
        # Silica fitting states
        self.silicaCA_fittingLine_drawn = False
        self.silica_autofit_done = False
        
        # Solvent CA fitting states
        self.solventCA_fittingLine_drawn = False
        self.solventCA_autofit_done = False
        
        # Solvent OA fitting states
        self.solventOA_fittingLine_drawn = False
        self.solventOA_fittingDone = False
        
        # Sample CA fitting states
        self.sampleCA_fittingLine_drawn = False
        self.sampleCA_autofit_done = False
        
        # Sample OA fitting states
        self.sampleOA_fittingLine_drawn = False
        self.sampleOA_autofit_done = False
        
        # Data processing states
        self.data_reversed = False
        self.where_to_start = None
    
    def reset_measurement_states(self):
        """Reset measurement-related states."""
        self.data_acquisition_complete = False
        self.data_reversed = False
        self.experiment_stopped = False
        self.running = False
    
    def reset_fitting_states(self, sample_type: str | None = None):
        """Reset fitting-related states.
        
        Args:
            sample_type: "silica", "solvent", or "sample". 
                        If None, resets all fitting states.
        """
        if sample_type is None or sample_type == "silica":
            self.silicaCA_fittingLine_drawn = False
            self.silica_autofit_done = False
        
        if sample_type is None or sample_type == "solvent":
            self.solventCA_fittingLine_drawn = False
            self.solventCA_autofit_done = False
            self.solventOA_fittingLine_drawn = False
            self.solventOA_fittingDone = False
        
        if sample_type is None or sample_type == "sample":
            self.sampleCA_fittingLine_drawn = False
            self.sampleCA_autofit_done = False
            self.sampleOA_fittingLine_drawn = False
            self.sampleOA_autofit_done = False
    
    def get_state_summary(self) -> dict:
        """Get a summary of all current states.
        
        Returns:
            Dictionary with all state variables and their current values.
        """
        return {
            "clearing": self.clearing,
            "data_acquisition_complete": self.data_acquisition_complete,
            "experiment_stopped": self.experiment_stopped,
            "running": self.running,
            "initialized": self.initialized,
            "initializing": self.initializing,
            "silicaCA_fittingLine_drawn": self.silicaCA_fittingLine_drawn,
            "silica_autofit_done": self.silica_autofit_done,
            "solventCA_fittingLine_drawn": self.solventCA_fittingLine_drawn,
            "solventCA_autofit_done": self.solventCA_autofit_done,
            "solventOA_fittingLine_drawn": self.solventOA_fittingLine_drawn,
            "solventOA_fittingDone": self.solventOA_fittingDone,
            "sampleCA_fittingLine_drawn": self.sampleCA_fittingLine_drawn,
            "sampleCA_autofit_done": self.sampleCA_autofit_done,
            "sampleOA_fittingLine_drawn": self.sampleOA_fittingLine_drawn,
            "sampleOA_autofit_done": self.sampleOA_autofit_done,
            "data_reversed": self.data_reversed,
            "where_to_start": self.where_to_start,
        }
