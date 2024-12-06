#include "FieldFactory.h"
#include "MultiPhaseField.h"
#include "LEBcMultiPhaseField.h"

FieldPtr FieldFactory::make_field() {
	// the default box is the cubic one
	char field_type[512] = "multiphasefield";

	if(strncmp(field_type, "multiphasefield", 512)) return std::make_shared<MultiPhaseField>();
	else if(strncmp(field_type, "lebcmultiphasefield", 512)) return std::make_shared<LEBcMultiPhaseField>();
	else throw("Unsupported field", field_type);
}
