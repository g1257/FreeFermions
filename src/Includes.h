#include "Engine.h"
#include "ObservableLibrary.h"
#include "GeometryLibrary.h"

typedef double RealType;
typedef double FieldType;
typedef std::vector<bool> LevelsType;
typedef psimag::Matrix<FieldType> MatrixType;
typedef FreeFermions::Engine<RealType,FieldType,LevelsType> EngineType;
typedef EngineType::HilbertVectorType HilbertVectorType;
typedef EngineType::FreeOperatorType FreeOperatorType;
typedef FreeFermions::GeometryLibrary<MatrixType> GeometryLibraryType;
typedef FreeFermions::ObservableLibrary<EngineType> ObservableLibraryType;
