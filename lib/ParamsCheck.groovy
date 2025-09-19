class ParamsCheck {
    def params
    def definitions

    ParamsCheck(def params, def definitions) {
        this.params = params
        this.definitions = definitions
    }

    def validateAndBuild() {
        def validated = [:]
        def config = this
        // First, update definitions with any user-provided params
        config.params.each { paramName, value ->
            if (paramName == 'definitions') return
            if (config.definitions.containsKey(paramName)) {
                config.definitions[paramName].default_value = value
            }
        }
        // Now, validate and build output
        config.definitions.each { paramName, defs ->
            def value = defs.default_value
            // Error if required but not provided
            if (defs.required && (value == null || (value instanceof String && value.trim().isEmpty()))) {
                throw new RuntimeException("Parameter '${paramName}' is required but was not provided.")
            }
            // If still no value, continue to next param
            if (value == null) {
                return
            }
            // Trim string values
            if (value instanceof String) {
                value = value.trim()
            }
            // Validate against 'allow' regex if defined
            if (defs.allow) {
                def matched = defs.allow.any { it == '*' || value.matches(it) }
                if (!matched) {
                    throw new RuntimeException("Value for parameter '${paramName}' does not match allowed patterns: '${value}'.")
                }
            }
            // Convert value to correct type if possible, error if conversion fails
            if (defs.type == 'integer') {
                try {
                    value = value as int
                } catch (Exception e) {
                    throw new RuntimeException("Parameter '${paramName}' expected an integer but got '${value}'.")
                }
                if (defs.min != null && value < defs.min) {
                    throw new RuntimeException("Value for parameter '${paramName}' must be >= ${defs.min}, but was ${value}.")
                }
                if (defs.max != null && value > defs.max) {
                    throw new RuntimeException("Value for parameter '${paramName}' must be <= ${defs.max}, but was ${value}.")
                }
            } else if (defs.type == 'float') {
                try {
                    value = value as float
                } catch (Exception e) {
                    throw new RuntimeException("Parameter '${paramName}' expected a float but got '${value}'.")
                }
                if (defs.min != null && value < defs.min) {
                    throw new RuntimeException("Value for parameter '${paramName}' must be >= ${defs.min}, but was ${value}.")
                }
                if (defs.max != null && value > defs.max) {
                    throw new RuntimeException("Value for parameter '${paramName}' must be <= ${defs.max}, but was ${value}.")
                }
            } else if (defs.type == 'path') {
                def file = new File(value)
                if (!file.exists()) {
                    throw new RuntimeException("Path for parameter '${paramName}' does not exist: '${value}'.")
                }
                value = file.toString()
            }
            // Optionally add more conversions here
            validated[paramName] = value
        }
        return validated
    }
}
